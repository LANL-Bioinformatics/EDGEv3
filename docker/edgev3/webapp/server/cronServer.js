const express = require('express');
require('dotenv').config();
const cors = require('cors');
const mongoose = require('mongoose');
const cron = require('node-cron');
const logger = require('./utils/logger');
const uploadMonitor = require('./crons/uploadMonitor');
const cromwellJobMonitor = require('./crons/cromwellJobMonitor');
const cromwellWorkflowMonitor = require('./crons/cromwellWorkflowMonitor');
const nextflowJobMonitor = require('./crons/nextflowJobMonitor');
const nextflowWorkflowMonitor = require('./crons/nextflowWorkflowMonitor');
const projectDeletionMonitor = require('./crons/projectDeletionMonitor');
const projectStatusMonitor = require('./crons/projectStatusMonitor');
const dbBackup = require('./crons/dbBackup');
const dbBackupClean = require('./crons/dbBackupClean');
const config = require('./config');

const app = express();
app.use(express.json());

// allow cross-origin requests
app.use(cors());

// cron jobs
// monitor cromwell jobs on every 2 minutes
cron.schedule(config.CRON.SCHEDULES.CROMWELL_JOB_MONITOR, async () => {
  await cromwellJobMonitor();
});
// monitor workflow requests on every 2 minutes
cron.schedule(config.CRON.SCHEDULES.CROMWELL_WORKFLOW_MONITOR, async () => {
  await cromwellWorkflowMonitor();
});
// monitor nextflow jobs on every 2 minutes
cron.schedule(config.CRON.SCHEDULES.NEXTFLOW_JOB_MONITOR, async () => {
  await nextflowJobMonitor();
});
// cron jobs
// monitor workflow requests on every 2 minutes
cron.schedule(config.CRON.SCHEDULES.NEXTFLOW_WORKFLOW_MONITOR, async () => {
  await nextflowWorkflowMonitor();
});
// monitor uploads every day at midnight
cron.schedule(config.CRON.SCHEDULES.FILE_UPLOAD_MONITOR, async () => {
  await uploadMonitor();
});
// monitor project status on every 1 minute
cron.schedule(config.CRON.SCHEDULES.PROJECT_STATUS_MONITOR, async () => {
  await projectStatusMonitor();
});
// monitor project deletion every day at 10pm
cron.schedule(config.CRON.SCHEDULES.PROJECT_DELETION_MONITOR, async () => {
  await projectDeletionMonitor();
});
// backup nmdcedge DB every day at 10pm
cron.schedule(config.CRON.SCHEDULES.DATABASE_BACKUP_CREATOR, () => {
  dbBackup();
});
// delete older DB backups every day at 12am
cron.schedule(config.CRON.SCHEDULES.DATABASE_BACKUP_PRUNER, () => {
  dbBackupClean();
});

const runApp = async () => {
  try {
    // Connect to MongoDB
    const db = `mongodb://${config.DATABASE.SERVER_HOST}:${config.DATABASE.SERVER_PORT}/${config.DATABASE.NAME}`;
    mongoose.set('strictQuery', false);
    mongoose
      .connect(
        db
      );
    logger.info(`Successfully connected to database ${db}`);
    // start server
    app.listen(config.CRON.SERVER_PORT, () => logger.info(`HTTP CRON server up and running on port ${config.CRON.SERVER_PORT} !`));
  } catch (err) {
    logger.error(err);
  }
};

runApp();
