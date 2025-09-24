const m2s = require('mongoose-to-swagger');
const User = require('../models/user');
const Project = require('../models/project');
const Upload = require('../models/upload');
const Job = require('../models/job');

module.exports = {
  user: { ...m2s(User) },
  /* registerUser: {
        ...m2s(User, {
            omitFields: ['_id', 'created', 'updated', 'job', 'notification', 'code', 'role'],
        }),
        ...registerUserExample,
    }, */
  project: m2s(Project, { omitFields: ['_id', 'created', 'updated', 'jobPriority', 'sharedTo', 'parent', 'children'] }),
  upload: m2s(Upload, { omitFields: ['_id', 'created', 'updated', 'sharedTo'] }),
  job: m2s(Job, { omitFields: ['_id', 'created', 'updated'] }),
};
