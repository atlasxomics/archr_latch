from wf import *

if __name__ == '__main__':
    archr_workflow(
        runs=[
            Run(
            'D01281_NG02604',
            LatchFile('latch://13502.account/cleaned/Wang_NCI/cleaned_D01281_NG02604_fragments.tsv.gz'),
            'control',
            LatchDir('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1281/spatial/'),
            LatchFile('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1281/spatial/tissue_positions_list.csv')
            ),
            Run(
            'D01282_NG02605',
            LatchFile('latch://13502.account/cleaned/Wang_NCI/cleaned_D01282_NG02605_fragments.tsv.gz'),
            'control',
            LatchDir('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1282/spatial/'),
            LatchFile('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1282/spatial/tissue_positions_list.csv')
            ),
            Run(
            'D01283_NG02608',
            LatchFile('latch://13502.account/cleaned/Wang_NCI/cleaned_D01283_NG02608_fragments.tsv.gz'),
            'control',
            LatchDir('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1283/spatial/'),
            LatchFile('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1283/spatial/tissue_positions_list.csv')
            ),
            Run(
            'D01284_NG12609',
            LatchFile('latch://13502.account/cleaned/Wang_NCI/cleaned_D01284_NG02609_fragments.tsv.gz'),
            'control',
            LatchDir('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1284/spatial/'),
            LatchFile('latch://13502.account/atx-illumina-1682977469.0200825/Images_spatial/D1284/spatial/tissue_positions_list.csv')
            ),
        ],
        project_name='multi_slims_test',
        genome=Genome.mm10,
        upload=True,
        lsi_resolution=[0.5],
        lsi_varfeatures=[25000],
        clustering_resolution=[1.0],
        min_TSS=1.5,
        min_frags=2500, 
    )
