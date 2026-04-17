package cz.siret.prank.program.routines.predict.external

import cz.siret.prank.program.PrankException
import cz.siret.prank.program.params.Parametrized
import cz.siret.prank.utils.Futils
import cz.siret.prank.utils.ProcessRunner
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Low-level runner for AutoDock Vina (or any Vina-compatible binary).
 *
 * Follows the same pattern as {@link FpocketRunner}: assembles the command line,
 * manages a per-run temp directory, redirects output to a log file, and raises a
 * {@link PrankException} on non-zero exit codes.
 */
@Slf4j
@CompileStatic
class VinaRunner implements Parametrized {

    /** Vina binary command (may include pre-set extra flags). */
    String vinaCommand = params.vina_command

    /**
     * Docking box specification.
     */
    static class BoxSpec {
        double centerX
        double centerY
        double centerZ
        double sizeX
        double sizeY
        double sizeZ

        BoxSpec(double cx, double cy, double cz, double sx, double sy, double sz) {
            this.centerX = cx; this.centerY = cy; this.centerZ = cz
            this.sizeX   = sx; this.sizeY   = sy; this.sizeZ   = sz
        }
    }

    /**
     * Result of a single Vina docking run.
     */
    static class VinaResult {
        /** Absolute path to the Vina output PDBQT file with all poses. */
        String outputPdbqt
        /** Absolute path to the Vina log/score text file. */
        String logFile
        /**
         * Parsed binding affinities (kcal/mol) for each pose, in pose order.
         * Empty if parsing failed.
         */
        List<Double> affinities = []
    }

    // -----------------------------------------------------------------------

    /**
     * Verify that the configured Vina binary is accessible and responds to --version.
     *
     * @throws PrankException if Vina is not found or returns a non-zero exit code.
     */
    void checkAvailable() {
        String cmd = "$vinaCommand --version"
        log.info "Checking Vina availability: [$cmd]"
        try {
            int code = ProcessRunner.process(cmd).redirectErrorStream().executeAndWait()
            if (code != 0) {
                throw new PrankException(
                    "Vina binary [$vinaCommand] returned exit code $code when invoked with --version. " +
                    "Make sure Vina is installed and the 'vina_command' parameter is correct.")
            }
        } catch (IOException e) {
            throw new PrankException(
                "Vina binary [$vinaCommand] could not be started: ${e.message}. " +
                "Make sure Vina is installed and available on PATH (or set 'vina_command' to the full path).", e)
        }
    }

    /**
     * Run Vina for a single receptor / ligand / box combination.
     *
     * @param receptorPdbqt  absolute path to the prepared receptor .pdbqt file
     * @param ligandPdbqt    absolute path to the prepared ligand .pdbqt file
     * @param box            docking box specification
     * @param outPdbqt       absolute path where Vina should write docked poses (.pdbqt)
     * @param tmpDir         working directory for this run (log file will be written here)
     * @return               {@link VinaResult} with paths and parsed affinities
     */
    VinaResult execute(String receptorPdbqt, String ligandPdbqt, BoxSpec box,
                       String outPdbqt, String tmpDir) {

        Futils.mkdirs(tmpDir)
        String logFile = "$tmpDir/vina.log"

        String command = buildCommand(receptorPdbqt, ligandPdbqt, box, outPdbqt)

        log.info "Running Vina: [$command] in [$tmpDir]"

        ProcessRunner runner = ProcessRunner.process(command, tmpDir)
                .redirectErrorStream()
                .redirectOutput(logFile)

        int exitcode = runner.executeAndWait()

        if (exitcode != 0) {
            throw new PrankException(
                "Vina failed with exit code $exitcode. " +
                "Receptor: [$receptorPdbqt], Ligand: [$ligandPdbqt]. " +
                "See Vina log: [$logFile]")
        }

        log.info "Vina finished successfully. Output: [$outPdbqt]"

        VinaResult result = new VinaResult()
        result.outputPdbqt = outPdbqt
        result.logFile     = logFile
        result.affinities  = parseAffinities(logFile)

        return result
    }

    // -----------------------------------------------------------------------

    private String buildCommand(String receptor, String ligand, BoxSpec box, String out) {
        StringBuilder sb = new StringBuilder(vinaCommand)
        sb.append(" --receptor ").append(receptor)
        sb.append(" --ligand "  ).append(ligand)
        sb.append(" --center_x ").append(format(box.centerX))
        sb.append(" --center_y ").append(format(box.centerY))
        sb.append(" --center_z ").append(format(box.centerZ))
        sb.append(" --size_x "  ).append(format(box.sizeX))
        sb.append(" --size_y "  ).append(format(box.sizeY))
        sb.append(" --size_z "  ).append(format(box.sizeZ))
        sb.append(" --out "     ).append(out)
        sb.append(" --num_modes ").append(params.vina_num_modes)
        sb.append(" --energy_range ").append(format(params.vina_energy_range))
        sb.append(" --exhaustiveness ").append(params.vina_exhaustiveness)
        if (params.vina_cpu > 0) {
            sb.append(" --cpu ").append(params.vina_cpu)
        }
        return sb.toString()
    }

    private static String format(double v) {
        String.format(Locale.US, "%.3f", v)
    }

    /**
     * Parse per-pose binding affinities from a Vina log file.
     *
     * Vina prints a table like:
     * <pre>
     *    mode |   affinity | dist from best mode
     *         | (kcal/mol) | rmsd l.b.| rmsd u.b.
     *    -----+------------+----------+----------
     *       1 |      -7.1  |   0.000  |   0.000
     *       2 |      -6.8  |   1.234  |   2.345
     * </pre>
     */
    static List<Double> parseAffinities(String logFile) {
        List<Double> affinities = []
        if (!Futils.exists(logFile)) return affinities

        boolean inTable = false
        try {
            new File(logFile).eachLine { String line ->
                if (line.contains('-----+------------')) {
                    inTable = true
                    return
                }
                if (inTable) {
                    // table lines: "   1 |      -7.1  |  0.000 | 0.000"
                    String trimmed = line.trim()
                    if (trimmed.empty || !Character.isDigit(trimmed.charAt(0))) {
                        inTable = false
                        return
                    }
                    String[] parts = trimmed.split(/\s*\|\s*/)
                    if (parts.length >= 2) {
                        try {
                            affinities << Double.parseDouble(parts[1].trim())
                        } catch (NumberFormatException ignored) {
                            // skip unparseable line
                        }
                    }
                }
            }
        } catch (Exception e) {
            log.warn "Could not parse Vina affinities from [$logFile]: ${e.message}"
        }
        return affinities
    }

}
