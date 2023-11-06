'''
Workflow for testing multiple ArchR parameter sets on ATAC fragments
files (fragments.tsv.gz); outputs umap plots, spatialdim plots, and
tss/fragment heatmaps.
'''

import glob
import re
import subprocess

from enum import Enum
from pathlib import Path
from typing import List
from latch.registry.table import Table

from latch import large_task, small_task, workflow, custom_task
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchParameter,
    LatchRule
)

import wf.lims as lims
from wf.registry import Run, upload_to_registry, Project, initialize_runs


class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'


@custom_task(cpu=62, memory=384, storage_gib=500)
# @large_task
def archr_task(
    projects: List[Project],
    project_name: str,
    genome: Genome,
    tile_size: int,
    min_TSS: float,
    min_frags: int,
    lsi_iterations: int,
    lsi_resolution: List[float],
    lsi_varfeatures: List[int],
    clustering_resolution: List[float],
    umap_mindist: float,
    project_table_id: str,
    run_table_id: str,
) -> LatchDir:

    runs =[]
    runs=initialize_runs(projects, project_table_id, run_table_id)
    _archr_cmd = [
        'Rscript',
        '/root/wf/archr_objs.R',
        project_name,
        genome.value,
        f'{tile_size}',
        f'{min_TSS}',
        f'{min_frags}',
        f'{lsi_iterations}',
        f'{",".join(str(i) for i in lsi_resolution)}',
        f'{",".join(str(i) for i in lsi_varfeatures)}',
        f'{",".join(str(i) for i in clustering_resolution)}',
        f'{umap_mindist}',
    ]

    runs = [
        (
            f'{run.run_id},'
            f'{run.fragments_file},'
            f'{run.condition},'
            f'{run.positions_file},'
            f'{run.spatial_dir},'
            )
        for run in runs
    ]

    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd)
    out_dir = project_name
    mkdir_cmd = ['mkdir'] + [out_dir]
    subprocess.run(mkdir_cmd)
    
    figures = glob.glob('*_plots.pdf')

    _mv_cmd = ['mv'] + figures + ['medians.csv'] + [out_dir]

    subprocess.run(_mv_cmd)

    return LatchDir(
        f'/root/{out_dir}',
        f'latch:///optimize_outs/{out_dir}'
    )


@small_task(retries=0)
def lims_task(
    results_dir: LatchDir,
    upload: bool,
) -> LatchDir:

    if upload:

        csv_path = Path(results_dir.local_path + '/medians.csv').resolve()
        slims = lims.slims_init()
        ng_re = re.compile('NG[0-9]{5}')

        with open(csv_path, 'r') as f:

            f.readline()
            lines = [
                tuple(line.split(',')) for line in
                [line.rstrip() for line in f.readlines()]
                ]

            for line in lines:
                run_id, tss, nfrags = line

                if re.findall(ng_re, run_id):

                    ng_id = re.findall(ng_re, run_id)[0]
                    print(f'Uploading results for {ng_id}')

                    try:
                        pk = lims.get_pk(ng_id, slims)
                    except IndexError:
                        print(f'Invalid SLIMS ng_id: {ng_id}.')
                        continue

                    payload = {}
                    payload['rslt_fk_content'] = pk
                    payload['rslt_fk_test'] = 55
                    payload['rslt_value'] = 'upload'
                    payload['rslt_cf_medianTssScore'] = tss
                    payload['rslt_cf_medianNfrags'] = nfrags

                    lims.push_result(payload, slims)
                    print("Upload to SLIMS succeeded.")

                else:
                    print(f"No NG_ID found for run {run_id}; upload failed.")

        return results_dir

    return results_dir


metadata = LatchMetadata(
    display_name='optimize archr',
    author=LatchAuthor(
        name='AtlasXomics, Inc.',
        email='jamesm@atlasxomics.com',
        github='github.com/atlasxomics',
    ),
    repository='https://github.com/atlasxomics/archr_latch',
    license='MIT',
    parameters={
        'projects': LatchParameter(
            display_name='projects',
            description='List of runs to be analyzed; each run must contain a \
                         run_id and fragments.tsv file; optional: condition, \
                         tissue position file for filtering on/off tissue, \
                         spatial folder for SpatialDimPlot.',
            batch_table_column=True,
            samplesheet=True
        ),
        'project_name': LatchParameter(
            display_name='project name',
            description='Name of output directory in optimize_outs/',
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*",
                    message="project name cannot start with a '/'"
                )
            ]
        ),
        'genome': LatchParameter(
            display_name='genome',
            description='Reference genome to be used for geneAnnotation and \
                        genomeAnnotation',
            batch_table_column=True,
        ),
        'upload': LatchParameter(
            display_name='upload to SLIMS',
            description='Upload median TSS and nFrags for each run to SLIMS; \
                        run_id MUST contain a valid NG ID (ie. NG12345).',
            batch_table_column=True,
        ),
        'tile_size': LatchParameter(
            display_name='tile size',
            description='The size of the tiles used for binning counts in the \
                        TileMatrix.',
            batch_table_column=True,
            hidden=True
        ),
        'min_TSS': LatchParameter(
            display_name='minimum TSS',
            description='The minimum numeric transcription start site (TSS) \
                        enrichment score required for a cell to pass \
                        filtering.',
            batch_table_column=True,
            hidden=True
        ),
        'min_frags': LatchParameter(
            display_name='minimum fragments',
            description='The minimum number of mapped ATAC-seq fragments \
                        required per cell to pass filtering.',
            batch_table_column=True,
            hidden=True
        ),
        'lsi_iterations': LatchParameter(
            display_name='LSI iterations',
            description='iterations parameter from addIterativeLSI function.',
            batch_table_column=True,
            hidden=True
        ),
        'lsi_resolution': LatchParameter(
            display_name='LSI resolution',
            description='resolution parameter from \
                        addIterativeLSI/clusterParams function.',
            batch_table_column=True
        ),
        'lsi_varfeatures': LatchParameter(
            display_name='LSI varFeatures',
            description='varFeatures parameter from addIterativeLSI function; \
                        each will correspond to a umap.pdf, the last in the \
                        will be used to make the RDS object.',
            batch_table_column=True
        ),
        'clustering_resolution': LatchParameter(
            display_name='clustering resolution',
            description='resolution parameter from addClusters function.',
            batch_table_column=True
        ),
        'umap_mindist': LatchParameter(
            display_name='UMAP minimum distance',
            description='minDist parameter from addUMAP function.',
            batch_table_column=True,
            hidden=True
        ),
        'run_table_id': LatchParameter(
            display_name='Runs Table ID',
            description='The ID of the runs table in Registry. \
            The runs will be updated in Registry with its \
                corresponding condition, spatial directory, condition, and \
                location of the optimized output archR project.'
        ),
        'project_table_id': LatchParameter(
            display_name='Projects Table ID',
            description='The ID of the projects/SOW table in Registry.\
            The optimized ArchR project will be inserted into \
                the SOW table for the corresponding runs.'
        )
    },
    tags=[],
)


@workflow(metadata)
def archr_workflow(
    projects: List[Project],
    genome: Genome,
    project_name: str,
    run_table_id: str = "761",
    project_table_id: str = "917",
    upload: bool = False,
    tile_size: int = 5000,
    min_TSS: float = 2.0,
    min_frags: int = 0,
    lsi_iterations: int = 2,
    lsi_resolution: List[float] = [0.5],
    lsi_varfeatures: List[int] = [25000],
    clustering_resolution: List[float] = [1.0],
    umap_mindist: float = 0.0
) -> LatchDir:
    ''' Workflow for assessing spatial epigenomic data generated via DBiT-seq.

    # optimize archr

    **optimize archr** is a [latch.bio]() workflow for assessing spatial
    epigenomic data generated via [DBiT-seq]
    (https://www.nature.com/articles/s41586-022-05094-1). Provided fragments
    from a single-cell ATAC-seq preprocessing and alignment
    workflow and spatial information, **optimize archr** returns plots and
    summary statistics to inform further processing.

    The workflow utilizes
    [ArchR](https://www.archrproject.com/articles/Articles/tutorial.html)
    to generate QC parameters (TSS, fragments per cell) and perform
    dimensionality reduction/clustering, and
    [Seurat](https://satijalab.org/seurat/)
    to spatially align the data.  The workflow can take data from either a
    single tissue-sample analyzed via DBiT-seq or multiple tissue-samples; in
    ATX parlance, tissue-samples analyzed via DBIT-seq are termed 'Runs'. Multiple Runs are linked
    to a Project, which is the input into **optimize archr**. All
    Runs inputted to **optimize archr** through Project(s) are merged into a single ArchRProject for
    analysis. 

    ## Inputs
    All input files for **optimize archr** must be on the latch.bio
    [file system](https://wiki.latch.bio/wiki/data/overview).  Each run in the
    workflow takes the following parameters,

    * [fragments.tsv.gz file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments):
    A BED-like, tab-delimited file in which each row contains an ATAC-seq
    fragment

    * [tissue_positions_list.csv](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html):
    A comma-separated file in which each row contains a unique barcode, an
    indicator for whether the tixel is 'on-tissue' (1, 0), and a row/column
    index

    * [Spatial folder](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html):
    A directory containing tissue images and experiment metadata

    * Run ID: An identifier for the run

    * Condition (_optional_):  An experimental Condition descriptor (ie.
    'control', 'diseased')

    Individual runs are linked to a Project, which are all batched with the following global
    parameters,

    * Project Name: A name for the output folder

    * Genome: A reference genome to be used for alignment

    * Upload to SLIMS _(ATX-internal only)_: A T/F toggle for whether to push
      QC results to LIMS

    * LSI resolution: A **list** of decimal values used as input to the
    `clusterParams` parameter of the `addIterativeLSI` function in
    [ArchR](https://www.archrproject.com/reference/addIterativeLSI.html);

    * LSI varFeatures: A **list** of integers used as input to the `varFeatures`
      parameter of the `addIterativeLSI` function in
      [ArchR](https://www.archrproject.com/reference/addIterativeLSI.html);

    * clustering resolution: A **list** of decimal values used as input to the
    `resolution` parameter of the `addClusters` function in
    [ArchR](https://www.archrproject.com/reference/addClusters.html).

    > The Project also takes a series of single-value parameters that can be
    found under the 'Hidden Parameters' dropdown; these parameters are less
    commonly varied inputs to ArchR functions.

    ## Running the workflow

    The **optimize archr** workflow can be found in the
    [Workflows](https://wiki.latch.bio/workflows/overview) module in your
    latch.bio workspace. For access to an ATX-collaborator workspace, please
    contact your AtlasXomics Support Scientist or email
    support@atlasxomics.com. See
    [here](https://wiki.latch.bio/workflows/overview) for general
    instructions for running workflows in latch.bio.

    1. Navigate to the **optimize archr** workflow in the Workflows module in
    your latch.bio workspace.  Ensure you are on the 'Parameters' tab of the
    workflow.

    2. To add Projects to the workflow, select the '+ Import from Registry' icon.  Select the Project(s)
    that are linked to the Runs you want to process. Repeat for each Project you want to add to the workflow.

    3. Scroll to the bottom of the page and input values for global project
    parameters.

    4. Click the 'Hidden Parameters' button and change the global parameters as
      needed.

    5. Click the 'Launch Workflow' button on the bottom-right of the parameters
    page.  This will automatically navigate you to the Executions tab of the
    workflow.

    6. From the Executions tab, you can view the status of the launched
    workflow. Once the workflow has completed running, the status will change
    to 'Succeeded'; if the workflow has the status 'Failed', please contact an
    AtlasXomics Support Scientist.  You can click on the workflow execution to
    view a more granular workflow status and see output logs.

    7. Workflow outputs are loaded into the latch.bio
    [Data module](https://wiki.latch.bio/wiki/data/overview) in the
    `optimize_outs` directory.

    ## Outputs

    Outputs from **optimize archr** are loaded into latch.bio
    [Data module](https://wiki.latch.bio/wiki/data/overview) in the
    `optimize_outs` directory.

    The workflow outputs four files:

    * medians.csv: A common-separated file containing median TSS and median
    fragment count values for on-tissue tixels for each Run in the Project

    * qc_plots.pdf: A PDF containing heat-maps for the fragment count (log10)
    and TSS Enrichment Score for each 'on-tissue' tixel, overlaid on top of
    the tissue image

    * umap_plots.pdf: A PDF containing
    [UMAP plots](https://www.archrproject.com/articles/Articles/tutorial.html#visualizing-in-a-2d-umap-embedding)
    colored by Run ID and cluster assignment. **optimize archr** generates a
    two UMAP plots for each element in the Cartesian product of the **LSI
    resolution**, **LSI varFeatures**, and **clustering resolution**
    parameters.

    * spatialdim_plots.pdf: A PDF containing Seurat
    [Spatial Plots](https://satijalab.org/seurat/reference/spatialplot) with
    tixel cluster-assignment overlaid on top of the tissue image.  For each
    Run, the workflow generates a Spatial Plot for each element in the
    Cartesian product of the **LSI resolution**, **LSI varFeatures**, and
    **clustering resolution** parameters.

    ## Next Steps

    The quality of each Run can be evaluated with the TSS Enrichment and
    Fragment count medians and the qc plots.  If the quality is deemed
    sufficient for further analysis, the UMAP plots and Spatial Plots can
    be used to inform dimensionality reduction and clustering in further
    analysis.

    Analysis can be performed locally or in a latch.bio
    [Pod](https://wiki.latch.bio/wiki/pods/overview).  For access to
    ATX-specific Pods, please contact your AtlasXomics Support Scientist.

    Further analysis can also be performed in latch.bio with the **create
    ArchRProject** (returns ArchRProject with peak and motif calling) and
    **atlasShiny** (returns inputs for the ATX ATAC-seq R Shiny App).  For
    access to these workflows, please contact your AtlasXomics Support
    Scientist.


    ## Support
    Questions? Comments?  Contact support@atlasxomics.com or post in
    AtlasXomics
    [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).

    '''

    results_dir = archr_task(
        projects=projects,
        project_name=project_name,
        genome=genome,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
        lsi_iterations=lsi_iterations,
        lsi_resolution=lsi_resolution,
        lsi_varfeatures=lsi_varfeatures,
        clustering_resolution=clustering_resolution,
        umap_mindist=umap_mindist,
        project_table_id=project_table_id,
        run_table_id=run_table_id
    )

    upload_to_registry(
        projects=projects,
        archr_project=results_dir,
        run_table_id=run_table_id,
        project_table_id=project_table_id
    )

    return lims_task(results_dir=results_dir, upload=upload)


LaunchPlan(
    archr_workflow,
    'defaults',
    {
    'projects' : [
        Project(
            'demo_row_archr', False
            )
        ],
    'project_name' : 'demo',
    'genome' : Genome.hg38,
    'run_table_id': '761',
    'project_table_id': '917',
    },
)

# if __name__ == "__main__":
#     archr_task(
#         projects=[Project(
#             project_id="example_proj_opt",
#             cleaned_frag_file=False
#         )],
#         project_name="test",
#         genome=Genome["mm10"],
#         tile_size=5000,
#         min_TSS=2.0,
#         min_frags=0,
#         lsi_iterations=2,
#         lsi_resolution=[0.5],
#         lsi_varfeatures=[25000],
#         clustering_resolution=[1.0],
#         umap_mindist=0.0,
#         project_table_id="922"
#     )