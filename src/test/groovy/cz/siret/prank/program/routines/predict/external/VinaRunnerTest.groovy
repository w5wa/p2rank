package cz.siret.prank.program.routines.predict.external

import groovy.transform.CompileStatic
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir

import java.nio.file.Path

import static org.junit.jupiter.api.Assertions.*

/**
 * Unit tests for {@link VinaRunner}.
 */
@CompileStatic
class VinaRunnerTest {

    // -----------------------------------------------------------------------
    // parseAffinities
    // -----------------------------------------------------------------------

    @Test
    void testParseAffinitiesTypical(@TempDir Path tmp) {
        String log = tmp.resolve("vina.log").toString()

        new File(log).text = """\
Detected 8 CPUs
Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 42

    mode |   affinity | dist from best mode
         | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
       1 |      -7.1  |   0.000  |   0.000
       2 |      -6.8  |   1.234  |   2.345
       3 |      -6.5  |   2.000  |   3.100
Writing output ... done.
"""

        List<Double> affinities = VinaRunner.parseAffinities(log)
        assertEquals(3, affinities.size())
        assertEquals(-7.1d, affinities[0], 0.001)
        assertEquals(-6.8d, affinities[1], 0.001)
        assertEquals(-6.5d, affinities[2], 0.001)
    }

    @Test
    void testParseAffinitiesEmptyFile(@TempDir Path tmp) {
        String log = tmp.resolve("empty.log").toString()
        new File(log).text = ""

        List<Double> affinities = VinaRunner.parseAffinities(log)
        assertTrue(affinities.empty)
    }

    @Test
    void testParseAffinitiesNoTable(@TempDir Path tmp) {
        String log = tmp.resolve("notable.log").toString()
        new File(log).text = "Vina crashed unexpectedly\nNo table here\n"

        List<Double> affinities = VinaRunner.parseAffinities(log)
        assertTrue(affinities.empty)
    }

    @Test
    void testParseAffinitiesMissingFile() {
        List<Double> affinities = VinaRunner.parseAffinities("/nonexistent/path/vina.log")
        assertTrue(affinities.empty)
    }

    // -----------------------------------------------------------------------
    // BoxSpec
    // -----------------------------------------------------------------------

    @Test
    void testBoxSpecCreation() {
        VinaRunner.BoxSpec box = new VinaRunner.BoxSpec(1.0, 2.0, 3.0, 20.0, 20.0, 20.0)
        assertEquals(1.0d, box.centerX, 0.0001)
        assertEquals(2.0d, box.centerY, 0.0001)
        assertEquals(3.0d, box.centerZ, 0.0001)
        assertEquals(20.0d, box.sizeX, 0.0001)
        assertEquals(20.0d, box.sizeY, 0.0001)
        assertEquals(20.0d, box.sizeZ, 0.0001)
    }

}
