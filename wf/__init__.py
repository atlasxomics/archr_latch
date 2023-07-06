''' 
Workflow for testing multiple ArchR parameter sets on ATAC fragments
files (fragments.tsv.gz); outputs umap plots, spatialdim plots, and
tss/fragment heatmaps.
'''

import glob
import re
import subprocess

from dataclasses import dataclass
from dataclasses_json import dataclass_json
from enum import Enum
from pathlib import Path
from typing import List

from latch import medium_task, small_task, workflow
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

@dataclass_json
@dataclass
class Run:
    run_id: str
    fragments_file: LatchFile
    condition: str = 'None'
    spatial_dir: LatchDir = LatchDir(
        'latch:///spatials/demo/spatial/'
    )
    positions_file: LatchFile = LatchFile(
        'latch:///spatials/demo/spatial/tissue_positions_list.csv'
    )

class Genome(Enum):
    mm10 = 'mm10'
    hg38 = 'hg38'

@medium_task
def archr_task(
    runs: List[Run],
    project_name: str,
    genome: Genome,
    tile_size: int,
    min_TSS: float,
    min_frags: int,
    lsi_iterations: int,
    lsi_resolution: List[float],
    lsi_varfeatures: List[int],
    clustering_resolution: List[float],
    umap_mindist: float
) -> LatchDir:
    
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
        f'{run.fragments_file.local_path},'
        f'{run.condition},'
        f'{run.positions_file.local_path},'
        f'{run.spatial_dir.local_path},'
        )
    for run in runs
    ]
    
    _archr_cmd.extend(runs)
    subprocess.run(_archr_cmd)

    out_dir = project_name
    subprocess.run(['mkdir', out_dir])

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
    
        csv_path = Path(results_dir.local_path + f'/medians.csv').resolve()
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
        'runs': LatchParameter(
            display_name='runs',
            description='List of runs to be analyzed; each run must contain a \
                         run_id and fragments.tsv file; optional: condition, \
                         tissue position file for filtering on/off tissue, \
                         spatial folder for SpatialDimPlot.',
            batch_table_column=True, 
        ),
        'project_name' : LatchParameter(
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
                        enrichment score required for a cell to pass filtering.',
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
    },
    tags=[],
)

@workflow(metadata)
def archr_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    upload: bool=False,
    tile_size: int=5000,
    min_TSS: float=2.0,
    min_frags: int=0,
    lsi_iterations: int=2,
    lsi_resolution: List[float]=[0.5],
    lsi_varfeatures: List[int]=[25000],
    clustering_resolution: List[float]=[1.0],
    umap_mindist: float=0.0
) -> LatchDir:
    '''Workflow for converting fragment.tsv.gz files from to ArchRProjects.

    For data from DBiT-seq for spatially-resolved epigenomics.
    - See Deng, Y. et al 2022.
    '''

    results_dir = archr_task(
        runs=runs,
        project_name=project_name,
        genome=genome,
        tile_size=tile_size,
        min_TSS=min_TSS,
        min_frags=min_frags,
        lsi_iterations=lsi_iterations,
        lsi_resolution=lsi_resolution,
        lsi_varfeatures=lsi_varfeatures,
        clustering_resolution=clustering_resolution,
        umap_mindist=umap_mindist
    )

    return lims_task(results_dir=results_dir, upload=upload)

LaunchPlan(
    archr_workflow,
    'defaults',
    {
    'runs' : [
        Run(
            'default',
            LatchFile('latch://13502.account/atac_outs/demo/outs/demo_fragments.tsv.gz'),
            'demo',
            LatchDir('latch://13502.account/spatials/demo/spatial'),
            LatchFile('latch://13502.account/spatials/demo/spatial/tissue_positions_list.csv'),
            )
        ],
    'project_name' : 'demo',
    'genome' : Genome.hg38,
    'upload' : False
    },
)