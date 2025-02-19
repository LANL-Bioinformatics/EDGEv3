export const workflowList = {
  sra2fastq: {
    label: 'Download SRA Data',
    category: 'data',
    // img: '/docs/images/sra2fastq.png',
    // thumbnail: '/docs/images/sra2fastq-thumbnail.png',
    link: 'https://github.com/LANL-Bioinformatics/EDGE_workflows/tree/main/sra2fastq',
    // doclink: 'https://nmdc-workflow-documentation.readthedocs.io/en/latest/chapters/6_MetaT_index.html',
    info: 'This tool retrieves sequence project in FASTQ files from NCBI- SRA / EBI - ENA / DDBJ database. Input accession number supports studies(SRP*/ ERP * /DRP*), experiments (SRX*/ERX * /DRX*), samples(SRS * /ERS*/DRS *), runs(SRR * /ERR*/DRR *), or submissions (SRA * /ERA*/DRA *).',
  },
  runFaQCs: {
    label: 'Reads QC',
    category: 'metagenomics',
    img: '/docs/images/runFaQCs.png',
    thumbnail: '/docs/images/runFaQCs-thumbnail.png',
    video: '/docs/videos/runFaQCs.mp4',
    pdf: '/docs/help/runFaQCs.pdf',
    link: 'https://github.com/microbiomedata/runFaQCs',
    doclink:
      'https://nmdc-workflow-documentation.readthedocs.io/en/latest/chapters/1_RQC_index.html',
    info: 'This workflow is a replicate of the QA protocol implemented at JGI for Illumina reads and use the program “rqcfilter2” from BBTools(38:44) which implements them as a pipeline.',
  },
}
