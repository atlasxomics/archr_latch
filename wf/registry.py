import logging

from latch import small_task
from latch.types import LatchDir, LatchFile
from latch.registry.table import Table
from dataclasses import dataclass
from typing import List
from latch.functions.messages import message

logging.basicConfig(format="%(levelname)s - %(asctime)s - %(message)s")

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

@small_task(retries=0)
def upload_to_registry(
    runs: List[Run],
    archr_project: LatchDir,
    run_table_id: str = "761",
    project_table_id: str = "779"
):
    run_table = Table(run_table_id)
    project_table = Table(project_table_id)
    try:
        with run_table.update() as updater:
            for run in runs:
                message(
                    "info",
                    {
                        "title": f"Updating run {run.run_id} in registry table ID {run_table_id}",
                        "body": f"Run {run.run_id}, condition {run.condition}",
                    },
                )

                updater.upsert_record(
                    run.run_id,
                    condition=run.condition,
                    spatial_directory=run.spatial_dir,
                    positions_file=run.positions_file,
                    optimize_outs=archr_project
                )
        
        with project_table.update() as updater:
            for run in runs:
                updater.upsert_record(
                    run.run_id,
                    optimize_outs=archr_project
                )
    except Exception as err:
        print(f"Unexpected {err=}, {type(err)=}")
        return
    finally:
        return


if __name__ == "__main__":
    upload_to_registry(
        runs=[
            Run(
                'D1291',
                LatchFile('latch:///atac_outs/6bp_D01291_NG02620/outs/6bp_D01291_NG02620_fragments.tsv.gz'),
                'old',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1291/spatial'),
                LatchFile('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1291/spatial/tissue_positions_list.csv'),
            ),
            Run(
                'D1292',
                LatchFile('latch:///cleaned/Babayev_cleaned/cleaned_D01292_NG02621_fragments.tsv.gz'),
                'old',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1292/spatial'),
                LatchFile('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1292/spatial/tissue_positions_list.csv')
            ),
            Run(
                'D1293',
                LatchFile('latch:///atac_outs/6bp_D01293_NG02622/outs/6bp_D01293_NG02622_fragments.tsv.gz'),
                'young',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1293/spatial'),
                LatchFile('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1293/spatial/tissue_positions_list.csv')
            ),
            Run(
                'D1294',
                LatchFile('latch:///atac_outs/D01294_NG02624/outs/D01294_NG02624_fragments.tsv.gz'),
                'young',
                LatchDir('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1294/spatial'),
                LatchFile('latch:///atx-illumina-1682977469.0200825/Images_spatial/D1294/spatial/tissue_positions_list.csv')
            )
        ],
        archr_project=LatchDir("latch://13502.account/ArchRProjects/Babeyev"),
        run_table_id="761",
        project_table_id="779"
    )
