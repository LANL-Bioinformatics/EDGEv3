const mongoose = require('mongoose');

// Create Schema
const userSchema = new mongoose.Schema({
  firstName: {
    type: String,
    required: true,
    trim: true,
  },
  lastName: {
    type: String,
    required: true,
    trim: true,
  },
  email: {
    type: String,
    required: true,
    unique: true,
    description: 'Same as user id',
    trim: true,
  },
  password: {
    type: String,
    required: true,
    description: 'At least 8 characters',
    trim: true,
  },
  code: {
    type: String,
    required: true,
    unique: false,
  },
  role: {
    type: String,
    default: 'user',
    enum: ['user', 'admin'],
  },
  active: {
    type: Boolean,
    default: false,
  },
  oauth: {
    type: String,
    required: false,
    description: 'Third-Party Authentication (orcid, google, facebook ...)',
    trim: true,
  },
  notification: {
    isOn: {
      type: Boolean,
      default: false,
    },
    email: {
      type: String,
      required: false,
      trim: true,
    },
  },
  job: {
    limit: {
      type: Number,
      default: 100,
    },
    priority: {
      type: Number,
      default: 0,
    },
  },
}, {
  timestamps: {
    createdAt: 'created', // Use `created` to store the created date
    updatedAt: 'updated', // and `updated` to store the last updated date
  },
}, {
  virtuals: {
    fullName: {
      get() {
        return `${this.firstName} ${this.lastName}`;
      },
    },
  },
});

module.exports = mongoose.model('User', userSchema);
