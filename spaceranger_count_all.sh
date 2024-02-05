########Spatial transcriptomic profiling of human retinoblastoma############
########created by Luozixian Wang and Raymond Wong##########################
########CERA, UNIMELB, 05/02/2024###########################################


## spaceranger count
# localcore # max 16 cores
# localmem # 128gb x 90% = 115gb

spaceranger count --id=RB_retina1 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20211203_slide1/mkfastq/HTHN5DRXY \
--sample=4_4m_A \
--image=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-A1_RB-1_C1533.jpg \
--slide=V10Y04-077  \
--area=A1 \
--loupe-alignment=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-A1.json \
--localcores=16 \
--localmem=115


spaceranger count --id=RB_retina1 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20211203_slide1/mkfastq/HTHN5DRXY \
--sample=4_4m_B \
--image=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-B1_D210190.jpg \
--slide=V10Y04-077  \
--area=B1 \
--loupe-alignment=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-B1.json \
--localcores=16 \
--localmem=115


spaceranger count --id=RB_retina1 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20211203_slide1/mkfastq/HTHN5DRXY \
--sample=4_4m_C \
--image=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-C1_D210305.jpg \
--slide=V10Y04-077  \
--area=C1 \
--loupe-alignment=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-C1.json \
--localcores=16 \
--localmem=115

spaceranger count --id=RB_retina1 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20211203_slide1/mkfastq/HTHN5DRXY \
--sample=4_4m_D \
--image=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-D1_D210028.jpg \
--slide=V10Y04-077  \
--area=D1 \
--loupe-alignment=/mnt/d/Visium_Rb/20211203_slide1/image/V10Y04-077-D1.json \
--localcores=16 \
--localmem=115

spaceranger count --id=V2A1_D210305S2 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20220228_slide2/mkfastq/HVWJMDRXY \
--sample=V2_4_4m_A1 \
--image=/mnt/d/Visium_Rb/20220228_slide2/image_stitch/V10Y04-078-A1.jpg \
--slide=V10Y04-078  \
--area=A1 \
--loupe-alignment=/mnt/d/Visium_Rb/20220228_slide2/image_stitch/V10Y04-078-A1.json \
--localcores=16 \
--localmem=115

spaceranger count --id=V2B1_RB_C1533_16 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20220228_slide2/mkfastq/HVWJMDRXY \
--sample=V2_4_4m_B1 \
--image=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-B1.jpg \
--slide=V10Y04-078  \
--area=B1 \
--loupe-alignment=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-B1.json \
--localcores=16 \
--localmem=115

spaceranger count --id=V2C1_RB_C1533_14 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20220228_slide2/mkfastq/HVWJMDRXY \
--sample=V2_4_4m_C1 \
--image=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-C1.jpg \
--slide=V10Y04-078  \
--area=C1 \
--loupe-alignment=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-C1.json \
--localcores=16 \
--localmem=115

spaceranger count --id=V2D1_RB_C1533_10 \
--transcriptome=/mnt/d/20210921_retina_snRNAseq_SP_S1/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/mnt/d/Visium_Rb/20220228_slide2/mkfastq/HVWJMDRXY \
--sample=V2_4_4m_D1 \
--image=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-D1.jpg \
--slide=V10Y04-078  \
--area=D1 \
--loupe-alignment=/mnt/d/Visium_Rb/20220228_slide2/image/V10Y04-078-D1.json \
--localcores=16 \
--localmem=115
