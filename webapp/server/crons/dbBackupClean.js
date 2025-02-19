const findRemoveSync = require('find-remove');

const logger = require('../utils/logger');
const config = require('../config');

module.exports = function dbBackup() {
  logger.debug('Clean up DB backup');
  const result = findRemoveSync(config.DATABASE.BACKUP_DIR, { dir: '^db-backup_', regex: true, age: { seconds: config.DATABASE.BACKUP_LIFETIME_SECONDS } });
  logger.info(result);
};
