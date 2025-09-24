// Note: first item will be the default value in DB
const projectStatus = ['in queue', 'running', 'failed', 'delete', 'rerun', 'interrupted', 'complete', 'processing', 'submitted'];
const uploadStatus = ['live', 'delete'];
const jobStatus = ['Submitted', 'Running', 'Failed', 'Aborted', 'Succeeded'];
const queueTypes = ['cromwell', 'nextflow'];

module.exports = {
  projectStatus,
  uploadStatus,
  jobStatus,
  queueTypes,
};
