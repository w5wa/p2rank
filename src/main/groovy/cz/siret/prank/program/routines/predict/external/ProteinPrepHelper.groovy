package cz.siret.prank.program.routines.predict.external

import cz.siret.prank.program.PrankException
import cz.siret.prank.program.params.Parametrized
import cz.siret.prank.utils.Futils
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Pre-docking preparation utilities.
 *
 * Cleans a PDB structure (removes waters, heteroatoms) and writes a minimal
 * Vina-readable PDBQT receptor by converting atom types.  This intentionally
 * avoids depending on external Python tools (OpenBabel, MGLTools) so that the
 * pipeline works without them.  Users who need full charge-assignment and
 * torsion-tree handling should supply a pre-prepared .pdbqt receptor directly
 * via the {@code vina_ligand} / receptor paths.
 *
 * <p>Metadata about the preparation (which filters were applied, source file,
 * timestamps) are written to a {@code prep_metadata.txt} file in the output
 * directory for reproducibility.
 */
@Slf4j
@CompileStatic
class ProteinPrepHelper implements Parametrized {

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Prepare a receptor PDB file for Vina docking.
     *
     * <ol>
     *   <li>Optionally removes HETATM (heteroatom) records.</li>
     *   <li>Optionally removes water molecules (HOH / WAT residues).</li>
     *   <li>Writes a cleaned .pdb file and a minimal .pdbqt copy.</li>
     *   <li>Writes {@code prep_metadata.txt} in outDir for reproducibility.</li>
     * </ol>
     *
     * @param sourcePdb  path to the input PDB file (may be compressed .gz)
     * @param outDir     directory where cleaned files will be written
     * @param label      short label used to construct output file names
     * @return           absolute path of the prepared .pdbqt receptor file
     */
    String prepareReceptor(String sourcePdb, String outDir, String label) {
        Futils.mkdirs(outDir)

        String cleanedPdb = "$outDir/${label}_cleaned.pdb"
        String receptorPdbqt = "$outDir/${label}_receptor.pdbqt"

        log.info "Preparing receptor: [$sourcePdb] -> [$receptorPdbqt]"

        // Read lines from (possibly compressed) source
        List<String> lines = readLines(sourcePdb)

        // Apply filters
        if (params.vina_prep_remove_water) {
            lines = removeWater(lines)
        }
        if (params.vina_prep_remove_heteroatoms) {
            lines = removeHeteroatoms(lines)
        }

        // Write cleaned PDB
        Futils.writeFile(cleanedPdb, lines.join('\n') + '\n')
        log.info "Cleaned PDB written: [$cleanedPdb]"

        // Convert to minimal PDBQT
        writePdbqt(lines, receptorPdbqt)
        log.info "Receptor PDBQT written: [$receptorPdbqt]"

        // Persist metadata
        writeMetadata(outDir, sourcePdb, cleanedPdb, receptorPdbqt)

        return Futils.absPath(receptorPdbqt)
    }

    /**
     * Validate that a ligand file exists and has a recognisable extension.
     * Raises {@link PrankException} on any problem.
     *
     * @param ligandPath path to the ligand file (.pdbqt, .sdf, .mol2)
     * @return           the absolute path (unchanged)
     */
    String validateLigand(String ligandPath) {
        if (!ligandPath || ligandPath.trim().empty) {
            throw new PrankException(
                "Ligand not specified. Set 'vina_ligand' to the path of a .pdbqt file.")
        }
        if (!Futils.exists(ligandPath)) {
            throw new PrankException("Ligand file not found: [$ligandPath]")
        }
        String ext = Futils.realExtension(ligandPath).toLowerCase()
        if (!(ext in ['pdbqt', 'sdf', 'mol2', 'pdb'])) {
            log.warn "Ligand file [$ligandPath] has extension '$ext' — " +
                     "Vina natively reads .pdbqt; other formats may require conversion."
        }
        return Futils.absPath(ligandPath)
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    private static List<String> readLines(String path) {
        if (Futils.isCompressed(path)) {
            return Futils.inputStream(path).text.readLines()
        }
        return new File(path).readLines()
    }

    private static List<String> removeWater(List<String> lines) {
        lines.findAll { String line ->
            // skip HETATM/ATOM records whose residue name is HOH or WAT
            if (!(line.startsWith('ATOM  ') || line.startsWith('HETATM'))) return true
            String resName = safeSubstring(line, 17, 20).trim().toUpperCase()
            return resName != 'HOH' && resName != 'WAT'
        }
    }

    private static List<String> removeHeteroatoms(List<String> lines) {
        lines.findAll { String line -> !line.startsWith('HETATM') }
    }

    /**
     * Write a minimal PDBQT file from cleaned PDB lines.
     *
     * AutoDock Vina can read PDBQT files without charge/torsion trees for the
     * rigid receptor; the key requirement is that columns 1-80 conform to PDB
     * format and that column 77-78 contains an AD4 atom type.  This method
     * appends a simple atom-type annotation (element symbol) to satisfy Vina.
     */
    private static void writePdbqt(List<String> pdbLines, String outPath) {
        StringBuilder sb = new StringBuilder()
        for (String line : pdbLines) {
            if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
                // PDB cols 77-78: AD4 atom type — use element symbol (cols 77-78 in 1-based)
                String padded = line.padRight(80)
                String element = safeSubstring(padded, 76, 78).trim()
                if (!element) {
                    // fallback: derive from atom name (col 13-14)
                    element = safeSubstring(padded, 12, 14).trim().replaceAll(/\d/, '').toUpperCase()
                    if (element.length() > 2) element = element[0..1]
                }
                // Build PDBQT line: original PDB record + partial charge (0.000) + atom type
                String pdbqtLine = padded[0..79] + String.format(Locale.US, " %6.3f", 0.0d) + " ${element}"
                sb.append(pdbqtLine).append('\n')
            } else if (line.startsWith('TER') || line.startsWith('END')) {
                sb.append(line).append('\n')
            } else if (!line.startsWith('REMARK') && !line.startsWith('HEADER') &&
                       !line.startsWith('TITLE')  && !line.startsWith('COMPND') &&
                       !line.startsWith('SOURCE') && !line.startsWith('AUTHOR') &&
                       !line.startsWith('REVDAT') && !line.startsWith('JRNL')   &&
                       !line.startsWith('DBREF')  && !line.startsWith('SEQRES') &&
                       !line.startsWith('SHEET')  && !line.startsWith('HELIX')) {
                // keep other records (CONECT, CRYST1 …) verbatim
                sb.append(line).append('\n')
            }
        }
        sb.append('END\n')
        Futils.writeFile(outPath, sb.toString())
    }

    private void writeMetadata(String outDir, String source, String cleanedPdb, String pdbqt) {
        String meta = [
            "source_pdb: $source",
            "cleaned_pdb: $cleanedPdb",
            "receptor_pdbqt: $pdbqt",
            "remove_water: ${params.vina_prep_remove_water}",
            "remove_heteroatoms: ${params.vina_prep_remove_heteroatoms}",
            "prepared_at: ${new Date()}",
        ].join('\n')
        Futils.writeFile("$outDir/prep_metadata.txt", meta + '\n')
    }

    private static String safeSubstring(String s, int start, int end) {
        if (s == null || s.length() <= start) return ''
        int safeEnd = Math.min(end, s.length())
        s.substring(start, safeEnd)
    }

}
