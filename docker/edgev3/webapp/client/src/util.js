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
    info: 'ReadsQC workflow ...',
  },
  assembly: {
    label: 'Assembly',
    category: 'metagenomics',
    info: 'Assembly workflow ...',
  },
  annotation: {
    label: 'Annotation',
    category: 'metagenomics',
    info: 'Annotation workflow ...',
  },
  binning: {
    label: 'Binning',
    category: 'metagenomics',
    info: 'Binning workflow ...',
  },
}
