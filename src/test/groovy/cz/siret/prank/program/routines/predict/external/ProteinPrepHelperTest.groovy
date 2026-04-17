package cz.siret.prank.program.routines.predict.external

import cz.siret.prank.program.PrankException
import cz.siret.prank.program.params.Params
import groovy.transform.CompileStatic
import org.junit.jupiter.api.BeforeEach
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir

import java.nio.file.Path

import static org.junit.jupiter.api.Assertions.*

/**
 * Unit tests for {@link ProteinPrepHelper}.
 */
@CompileStatic
class ProteinPrepHelperTest {

    // Minimal PDB content used as test input
    static final String MINIMAL_PDB = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
HETATM    3  O   HOH A 100       5.000   5.000   5.000  1.00  0.00           O
HETATM    4  FE  HEM A 200       6.000   6.000   6.000  1.00  0.00          FE
END
"""

    @BeforeEach
    void resetParams() {
        Params.INSTANCE = new Params()
    }

    // -----------------------------------------------------------------------
    // prepareReceptor – water removal
    // -----------------------------------------------------------------------

    @Test
    void testPrepareReceptorRemovesWater(@TempDir Path tmp) {
        File pdb = tmp.resolve("protein.pdb").toFile()
        pdb.text = MINIMAL_PDB

        Params.inst.vina_prep_remove_water = true
        Params.inst.vina_prep_remove_heteroatoms = false

        ProteinPrepHelper helper = new ProteinPrepHelper()
        String pdbqt = helper.prepareReceptor(pdb.absolutePath, tmp.resolve("prep").toString(), "test")

        assertTrue(new File(pdbqt).exists(), "PDBQT file should be created")
        String content = new File(pdbqt).text
        assertFalse(content.contains('HOH'), "Water should have been removed")
    }

    // -----------------------------------------------------------------------
    // prepareReceptor – heteroatom removal
    // -----------------------------------------------------------------------

    @Test
    void testPrepareReceptorRemovesHeteroatoms(@TempDir Path tmp) {
        File pdb = tmp.resolve("protein.pdb").toFile()
        pdb.text = MINIMAL_PDB

        Params.inst.vina_prep_remove_water = false
        Params.inst.vina_prep_remove_heteroatoms = true

        ProteinPrepHelper helper = new ProteinPrepHelper()
        String pdbqt = helper.prepareReceptor(pdb.absolutePath, tmp.resolve("prep").toString(), "test")

        String content = new File(pdbqt).text
        assertFalse(content.contains('HETATM'), "HETATM records should have been removed")
        assertTrue(content.contains('ATOM  '), "ATOM records should be retained")
    }

    // -----------------------------------------------------------------------
    // prepareReceptor – metadata file
    // -----------------------------------------------------------------------

    @Test
    void testPrepareReceptorWritesMetadata(@TempDir Path tmp) {
        File pdb = tmp.resolve("protein.pdb").toFile()
        pdb.text = MINIMAL_PDB

        Params.inst.vina_prep_remove_water = true
        Params.inst.vina_prep_remove_heteroatoms = true

        String prepDir = tmp.resolve("prep").toString()
        ProteinPrepHelper helper = new ProteinPrepHelper()
        helper.prepareReceptor(pdb.absolutePath, prepDir, "test")

        File meta = new File("$prepDir/prep_metadata.txt")
        assertTrue(meta.exists(), "Metadata file should be created")
        assertTrue(meta.text.contains("remove_water: true"))
        assertTrue(meta.text.contains("remove_heteroatoms: true"))
    }

    // -----------------------------------------------------------------------
    // validateLigand
    // -----------------------------------------------------------------------

    @Test
    void testValidateLigandThrowsForMissingFile() {
        ProteinPrepHelper helper = new ProteinPrepHelper()
        assertThrows(PrankException.class) {
            helper.validateLigand("/nonexistent/ligand.pdbqt")
        }
    }

    @Test
    void testValidateLigandThrowsForBlankPath() {
        ProteinPrepHelper helper = new ProteinPrepHelper()
        assertThrows(PrankException.class) {
            helper.validateLigand("   ")
        }
    }

    @Test
    void testValidateLigandReturnsAbsPath(@TempDir Path tmp) {
        File lig = tmp.resolve("ligand.pdbqt").toFile()
        lig.text = "REMARK ligand\n"
        ProteinPrepHelper helper = new ProteinPrepHelper()
        String result = helper.validateLigand(lig.absolutePath)
        assertEquals(lig.absolutePath, result)
    }

}
