package cz.siret.prank.program.routines.predict

import cz.siret.prank.domain.Dataset
import cz.siret.prank.domain.Pocket
import cz.siret.prank.domain.Prediction
import cz.siret.prank.domain.PredictionPair
import cz.siret.prank.features.FeatureExtractor
import cz.siret.prank.prediction.pockets.rescorers.ModelBasedRescorer
import cz.siret.prank.program.PrankException
import cz.siret.prank.program.ml.Model
import cz.siret.prank.program.routines.Routine
import cz.siret.prank.program.routines.predict.external.ProteinPrepHelper
import cz.siret.prank.program.routines.predict.external.VinaRunner
import cz.siret.prank.program.routines.predict.external.VinaRunner.BoxSpec
import cz.siret.prank.program.routines.predict.external.VinaRunner.VinaResult
import cz.siret.prank.utils.Futils
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.biojava.nbio.structure.Atom

import static cz.siret.prank.utils.ATimer.startTimer
import static cz.siret.prank.utils.Futils.mkdirs
import static cz.siret.prank.utils.Futils.writeFile

/**
 * Routine for the 'vina-dock' command.
 *
 * <p>Supports two modes controlled by {@code vina_use_p2rank_pockets}:</p>
 * <ul>
 *   <li><b>P2Rank-guided</b> (default): runs P2Rank pocket prediction first, then
 *       docks the ligand into the top-N predicted binding-site boxes.</li>
 *   <li><b>Normal Vina</b>: uses the explicit box defined by
 *       {@code vina_center_x/y/z} and {@code vina_size_x/y/z}.</li>
 * </ul>
 *
 * <p>For each protein in the dataset the routine:</p>
 * <ol>
 *   <li>Prepares the receptor (.pdbqt) via {@link ProteinPrepHelper}.</li>
 *   <li>Validates the ligand file.</li>
 *   <li>Runs Vina (once per box / pocket).</li>
 *   <li>Writes per-protein and per-pocket CSV summaries and keeps/deletes
 *       intermediate files based on {@code vina_keep_output}.</li>
 * </ol>
 */
@Slf4j
@CompileStatic
class VinaDockingRoutine extends Routine {

    Dataset dataset
    String modelf   // P2Rank model file (used only in P2Rank-guided mode)

    VinaDockingRoutine(Dataset dataset, String modelf, String outdir) {
        super(outdir)
        this.dataset = dataset
        this.modelf = modelf
    }

    // -----------------------------------------------------------------------

    Dataset.Result execute() {
        def timer = startTimer()

        mkdirs(outdir)
        writeParams(outdir)

        log.info "outdir: $outdir"

        write "Starting vina-dock for proteins from dataset [${dataset.name}]"
        write "Mode: ${params.vina_use_p2rank_pockets ? 'P2Rank-guided' : 'Normal Vina'}"

        // Validate Vina availability early
        new VinaRunner().checkAvailable()

        // Validate ligand early
        String ligandFile = new ProteinPrepHelper().validateLigand(params.vina_ligand)
        write "Ligand: [$ligandFile]"

        // Load P2Rank model only in guided mode
        Model model = null
        FeatureExtractor extractor = null
        if (params.vina_use_p2rank_pockets) {
            model    = Model.load(modelf)
            extractor = FeatureExtractor.createFactory()
        }

        final Model finalModel     = model
        final FeatureExtractor finalExtractor = extractor

        Dataset.Result result = dataset.processItems { Dataset.Item item ->

            String itemLabel = item.label
            String itemOutDir = "$outdir/$itemLabel"
            mkdirs(itemOutDir)

            write "\n--- Docking protein [$itemLabel] ---"

            // 1. Prepare receptor
            ProteinPrepHelper prepHelper = new ProteinPrepHelper()
            String receptorPdbqt = prepHelper.prepareReceptor(item.proteinFile, "$itemOutDir/prep", itemLabel)

            // 2. Determine box(es)
            List<BoxSpec> boxes
            List<String>  pocketLabels

            if (params.vina_use_p2rank_pockets) {
                // P2Rank-guided: predict pockets and use their centroids as box centers
                List predictResult = predictBoxes(item, finalModel, finalExtractor)
                boxes        = (List<BoxSpec>)  predictResult[0]
                pocketLabels = (List<String>)   predictResult[1]
                if (boxes.empty) {
                    throw new PrankException(
                        "P2Rank found no pockets in protein [$itemLabel]. " +
                        "Try increasing 'vina_max_pockets' or switch to normal Vina mode.")
                }
                write "Using ${boxes.size()} P2Rank pocket(s) as docking box center(s)"
            } else {
                boxes       = [new BoxSpec(
                    params.vina_center_x, params.vina_center_y, params.vina_center_z,
                    params.vina_size_x,   params.vina_size_y,   params.vina_size_z)]
                pocketLabels = ['box1']
                write "Using explicit box: center=(${params.vina_center_x},${params.vina_center_y},${params.vina_center_z})"
            }

            // 3. Run Vina for each box
            List<Map<String, Object>> allResults = []
            VinaRunner runner = new VinaRunner()

            boxes.eachWithIndex { BoxSpec box, int idx ->
                String pocketLabel = pocketLabels[idx]
                String runDir      = "$itemOutDir/docking_${pocketLabel}"
                mkdirs(runDir)

                String outPdbqt = "$runDir/docked_poses.pdbqt"

                write "  Running Vina for pocket [$pocketLabel] ..."

                VinaResult vinaResult = runner.execute(receptorPdbqt, ligandFile, box, outPdbqt, runDir)

                Map<String, Object> row = buildResultRow(itemLabel, pocketLabel, box, vinaResult)
                allResults << row

                write "  Pocket [$pocketLabel]: best affinity = ${vinaResult.affinities ? vinaResult.affinities[0] : 'N/A'} kcal/mol"
            }

            // 4. Write summary CSV
            String csvFile = "$itemOutDir/${itemLabel}_docking_results.csv"
            writeFile(csvFile, buildCsv(allResults))
            write "Results written: [$csvFile]"

            // 5. Clean up if requested
            if (!params.vina_keep_output) {
                log.info "Deleting intermediate docking files for [$itemLabel]"
                Futils.delete("$itemOutDir/prep")
            }
        }

        write "\nvina-dock finished in $timer.formatted"
        write "Results saved to [${Futils.absPath(outdir)}]"

        return result
    }

    // -----------------------------------------------------------------------
    // P2Rank pocket prediction
    // -----------------------------------------------------------------------

    /**
     * Run P2Rank on a single item and return box specs derived from pocket centroids.
     *
     * @return a two-element list: [List&lt;BoxSpec&gt;, List&lt;String&gt; pocketLabels]
     */
    private List predictBoxes(Dataset.Item item, Model model, FeatureExtractor extractor) {
        PredictionPair pair    = item.predictionPair
        ModelBasedRescorer rsc = new ModelBasedRescorer(model, extractor)
        rsc.reorderPockets(pair.prediction, item.context)

        Prediction prediction = pair.prediction
        List<Pocket> topPockets = prediction.pockets.take(params.vina_max_pockets)

        List<BoxSpec>  boxes  = []
        List<String>   labels = []

        topPockets.eachWithIndex { Pocket pocket, int idx ->
            Atom c = pocket.centroid
            if (c == null) {
                log.warn "Pocket ${idx + 1} has no centroid — skipping"
                return
            }
            boxes  << new BoxSpec(c.x, c.y, c.z,
                                  params.vina_size_x, params.vina_size_y, params.vina_size_z)
            labels << ("pocket${idx + 1}".toString())
        }

        return [boxes, labels]
    }

    // -----------------------------------------------------------------------
    // CSV output
    // -----------------------------------------------------------------------

    private static Map<String, Object> buildResultRow(String protein, String pocket,
                                                       BoxSpec box, VinaResult vr) {
        List<Double> affinities = vr.affinities
        [
            protein       : protein,
            pocket        : pocket,
            box_center_x  : box.centerX,
            box_center_y  : box.centerY,
            box_center_z  : box.centerZ,
            box_size_x    : box.sizeX,
            box_size_y    : box.sizeY,
            box_size_z    : box.sizeZ,
            num_poses     : affinities.size(),
            best_affinity : affinities ? affinities[0] : Double.NaN,
            all_affinities: affinities.join(';'),
            output_pdbqt  : vr.outputPdbqt,
            log_file      : vr.logFile,
        ] as Map<String, Object>
    }

    private static String buildCsv(List<Map<String, Object>> rows) {
        if (rows.empty) return ''

        List<String> headers = rows[0].keySet().toList()
        StringBuilder sb = new StringBuilder()
        sb.append(headers.join(',') + '\n')

        rows.each { Map<String, Object> row ->
            List<String> vals = headers.collect { String h ->
                def v = row.get(h)
                (v == null) ? '' : formatCsvValue(v.toString())
            }
            sb.append(vals.join(',') + '\n')
        }
        return sb.toString()
    }

    private static String formatCsvValue(String value) {
        if (value.contains(',') || value.contains('"') || value.contains('\n')) {
            '"' + value.replace('"', '""') + '"'
        } else {
            value
        }
    }

}
