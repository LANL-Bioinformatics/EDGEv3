const path = require('path');
const fs = require('fs');
const Papa = require('papaparse');
const Upload = require('../edge-api/models/upload');
const config = require('../config');

const cromwellWorkflows = [];
const nextflowWorkflows = [
  'sra2fastq',
  'runFaQCs',
  'assembly',
  'annotation',
  'binning',
  'antiSmash',
  'taxonomy',
  'phylogeny',
  'refBased',
  'geneFamily'
];
const nextflowConfigs = {
  executor_config: {
    slurm: 'slurm.config',
    local: 'local.config',
  },
  module_params: 'module_params.tmpl',
  container_config: 'container.config',
  nf_reports: 'nf_reports.tmpl',
};

const workflowList = {
  default_wdl_version: '1.0',
  sra2fastq: {
    // cromwell
    // set if not default 1.0
    // wdl_version: '1.0',
    wdl: 'data/sra2fastq.wdl',
    wdl_imports: 'data/imports.zip',
    inputs_tmpl: 'data/sra2fastq_inputs.tmpl',
    cromwell_calls: ['sra.sra2fastq'],
    outdir: 'output/sra2fastq',
    // nextflow
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',

  },
  runFaQCs: {
    outdir: 'output/ReadsQC',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  assembly: {
    outdir: 'output/Assembly',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  annotation: {
    outdir: 'output/Annotation',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  binning: {
    outdir: 'output/Binning',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  antiSmash: {
    outdir: 'output/AntiSmash',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  taxonomy: {
    outdir: 'output/Taxonomy',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  phylogeny: {
    outdir: 'output/Phylogeny',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  refBased: {
    outdir: 'output/RefBased',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
  geneFamily: {
    outdir: 'output/GeneFamily',
    nextflow_main: 'main.nf',
    config_tmpl: 'workflow_config.tmpl',
  },
};

const linkUpload = async (fq, projHome) => {
  try {
    if (fq.startsWith(config.IO.UPLOADED_FILES_DIR)) {
      // create input dir and link uploaded file with realname
      const inputDir = `${projHome}/input`;
      if (!fs.existsSync(inputDir)) {
        fs.mkdirSync(inputDir);
      }
      const fileCode = path.basename(fq);
      let name = fileCode;
      const upload = await Upload.findOne({ 'code': name });
      if (upload) {
        name = upload.name;
      }
      let linkFq = `${inputDir}/${name}`;
      let i = 1;
      while (fs.existsSync(linkFq)) {
        i += 1;
        if (name.includes('.')) {
          const newName = name.replace('.', `${i}.`);
          linkFq = `${inputDir}/${newName}`;
        } else {
          linkFq = `${inputDir}/${name}${i}`;
        }
      }
      fs.symlinkSync(fq, linkFq, 'file');
      return linkFq;
    }
    return fq;
  } catch (err) {
    return Promise.reject(err);
  }
};

const generateWorkflowResult = (proj) => {
  const projHome = `${config.IO.PROJECT_BASE_DIR}/${proj.code}`;
  const resultJson = `${projHome}/result.json`;

  if (!fs.existsSync(resultJson)) {
    const result = {};
    const projectConf = JSON.parse(fs.readFileSync(`${projHome}/conf.json`));
    const outdir = `${projHome}/${workflowList[projectConf.workflow.name].outdir}`;

    if (projectConf.workflow.name === 'sra2fastq') {
      // use relative path
      const { accessions } = projectConf.workflow.input;
      accessions.forEach((accession) => {
        // link sra downloads to project output
        fs.symlinkSync(`../../../../sra/${accession}`, `${outdir}/${accession}`);

      });
    } else if (projectConf.workflow.name === 'runFaQCs') {
      const statsJsonFile = `${outdir}/QC.stats.json`;
      if (fs.existsSync(statsJsonFile)) {
        result.stats = JSON.parse(fs.readFileSync(statsJsonFile));
      }
      const summaryPlotsFile = `${outdir}/QC_summary_plots.html`;
      if (fs.existsSync(summaryPlotsFile)) {
        result.summaryPlots = `${workflowList[projectConf.workflow.name].outdir}/QC_summary_plots.html`;
      }
      const reportFile = `${outdir}/QC_final_report.html`;
      if (fs.existsSync(reportFile)) {
        result.report = `${workflowList[projectConf.workflow.name].outdir}/QC_final_report.html`;
      }
      const reportLongReadsFile = `${outdir}/NanoPlot-report.html`;
      if (fs.existsSync(reportLongReadsFile)) {
        result.report = `${workflowList[projectConf.workflow.name].outdir}/NanoPlot-report.html`;
      }
    } else if (projectConf.workflow.name === 'assembly') {
      const statsFile = `${outdir}/contigs_stats.txt`;
      if (fs.existsSync(statsFile)) {
        result.stats = Papa.parse(fs.readFileSync(statsFile).toString(), { delimiter: '\t', header: true, skipEmptyLines: true }).data;
      }
      const reportFile = `${outdir}/final_report.pdf`;
      if (fs.existsSync(reportFile)) {
        result.report = `${workflowList[projectConf.workflow.name].outdir}/final_report.pdf`;
      }
    } else if (projectConf.workflow.name === 'phylogeny') {
      const treeAllHtml = `${outdir}/SNPphyloTree.all.html`;
      if (fs.existsSync(treeAllHtml)) {
        result.treeAllHtml = `${workflowList[projectConf.workflow.name].outdir}/SNPphyloTree.all.html`;
      }
      const treeCdsHtml = `${outdir}/SNPphyloTree.cds.html`;
      if (fs.existsSync(treeCdsHtml)) {
        result.treeCdsHtml = `${workflowList[projectConf.workflow.name].outdir}/SNPphyloTree.cds.html`;
      }
    } else if (projectConf.workflow.name === 'antiSmash') {
      // antiSMASH HTML output
      const antiSmashHtml = `${outdir}/output/index.html`;
      if (fs.existsSync(antiSmashHtml)) {
        result.antiSmashHtml = `${workflowList[projectConf.workflow.name].outdir}/output/index.html`;
      }
    }

    fs.writeFileSync(resultJson, JSON.stringify(result));
  }
};

module.exports = {
  cromwellWorkflows,
  nextflowWorkflows,
  nextflowConfigs,
  workflowList,
  linkUpload,
  generateWorkflowResult,
};
