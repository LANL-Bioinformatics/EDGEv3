const mongoose = require('mongoose');
const { projectStatus } = require('../utils/conf');

// Create Schema
const projectSchema = new mongoose.Schema({
  name: {
    type: String,
    required: true,
    trim: true,
  },
  desc: {
    type: String,
    trim: true,
  },
  status: {
    type: String,
    default: projectStatus[0],
    enum: projectStatus,
  },
  type: {
    type: String,
    required: true,
  },
  code: {
    type: String,
    required: true,
    unique: true,
  },
  public: {
    type: Boolean,
    default: false,
  },
  notified: {
    type: Boolean,
    default: false,
  },
  owner: {
    type: String,
    required: true,
  },
  jobPriority: {
    type: Number,
    default: 0,
  },
  sharedTo: [{
    type: String,
  }],
  parent: {
    type: String,
    required: false,
  },
  children: [{
    type: String,
  }],
}, {
  timestamps: {
    createdAt: 'created', // Use `created` to store the created date
    updatedAt: 'updated', // and `updated` to store the last updated date
  },
});

module.exports = mongoose.model('Project', projectSchema);
