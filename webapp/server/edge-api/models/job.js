const mongoose = require('mongoose');
const { jobStatus, queueTypes } = require('../utils/conf');

// Create Schema
const jobSchema = new mongoose.Schema({
  id: {
    type: String,
    required: true,
  },
  project: {
    type: String,
    required: true,
    unique: true,
  },
  type: {
    type: String,
    required: true,
  },
  inputSize: {
    type: Number,
    default: 0,
  },
  status: {
    type: String,
    default: jobStatus[0],
    enum: jobStatus,
  },
  queue: {
    type: String,
    default: queueTypes[0],
    enum: queueTypes,
  },
}, {
  timestamps: {
    createdAt: 'created', // Use `created` to store the created date
    updatedAt: 'updated', // and `updated` to store the last updated date
  },
});

module.exports = mongoose.model('Job', jobSchema);
