const registerUserExample = {
  example: {
    firstName: 'edge', lastName: 'lanl', email: 'my@edge.lanl', password: 'my#4EDGE', confirmPassword: 'my#4EDGE', active: true,
  },
};
const activateUserExample = { example: { email: 'my@edge.lanl', token: '$2a$10$5S/XP1bN0BOdj/MSYNhPO.Ee1TCb176IylCgDM5nz6dZ0vuKDr3Bu' } };
const resetPasswordExample = { example: { email: 'my@edge.lanl', token: '$2a$10$5S/XP1bN0BOdj/MSYNhPO.Ee1TCb176IylCgDM5nz6dZ0vuKDr3Bu', newPassword: 'my#4EDGE' } };
const getActionLinkExample = { example: { email: 'my@edge.lanl', actionURL: 'http://my.edge/activate' } };
const loginExample = { example: { email: 'my@edge.lanl', password: 'my#4EDGE' } };
const oauthLoginExample = {
  example: {
    firstName: 'edge', lastName: 'lanl', email: '0000-1111-1111-1111@orcid.org', oauth: 'orcid', token: 'cAjCpnNzbZQ4V3GU',
  },
};
const updateUserExample = {
  example: {
    firstName: 'edge', lastName: 'lanl', password: 'my#4EDGE', confirmPassword: 'my#4EDGE',
  },
};
const addProjectExample = {
  example: {
    project: {
      name: 'project-name', desc: 'project-desc',
    },
    pipeline: 'edge pipeline',
  },
};
const updateProjectExample = {
  example: {
    name: 'newproject-name', desc: 'newproject-desc', public: false, sharedTo: ['test@my.edge'],
  },
};
const addUploadExample = {
  example: {
    name: 'test.fa',
    desc: 'test',
    type: 'fa',
    size: 150
  },
};
const updateUploadExample = {
  example: {
    name: 'newUpload-name', desc: 'newUpload-desc', public: true, sharedTo: ['test@my.edge'],
  },
};
const projectFilesExample = {
  example: {
    fileTypes: ['fa', 'fa.gz', 'fasta', 'fasta.gz'],
  },
};
const projectActionExample = { example: { code: 'some' } };
const actionSuccessfulExample = { example: { message: 'Action successful', success: true } };
const actionFailedExample = { example: { message: 'Action failed', success: false } };
const serverErrorExample = { example: { message: 'API server error', success: false } };

module.exports = {
  registerUser: {
    properties: {
      firstName: {
        description: 'user first name', type: 'string',
      },
      lastName: {
        description: 'user last name', type: 'string',
      },
      email: {
        description: 'user email', type: 'string',
      },
      password: {
        description: 'password', type: 'string',
      },
      confirmPassword: {
        description: 'confirm password', type: 'string',
      },
      active: {
        description: 'user status, default false', type: 'boolean',
      },
      actionURL: {
        description: 'the URL of the activation page', type: 'string',
      },
    },
    required: ['firstName', 'lastName', 'email', 'password', 'confirmPassword'],
    ...registerUserExample,
  },
  activateUser: {
    properties: {
      email: {
        description: 'user email', type: 'string',
      },
      token: {
        description: 'encoded password', type: 'string',
      },
    },
    required: ['email', 'token'],
    ...activateUserExample,
  },
  resetPassword: {
    properties: {
      email: {
        description: 'user email', type: 'string',
      },
      token: {
        description: 'encoded password', type: 'string',
      },
      newPassword: {
        description: 'new password', type: 'string',
      },
    },
    required: ['email', 'token'],
    ...resetPasswordExample,
  },
  getActionLink: {
    properties: {
      email: {
        description: 'user email', type: 'string',
      },
      actionURL: {
        description: 'the URL of the client side action page', type: 'string',
      },
    },
    required: ['email'],
    ...getActionLinkExample,
  },
  login: {
    properties: {
      email: {
        description: 'user email', type: 'string',
      },
      password: {
        description: 'password', type: 'string',
      },
    },
    required: ['email', 'password'],
    ...loginExample,
  },
  oauthLogin: {
    properties: {
      email: {
        description: 'user email', type: 'string',
      },
      oauth: {
        description: 'Third-Party Authentication (orcid, google, facebook ...)', type: 'string',
      },
      token: {
        description: 'Our secret key. You have to obtain this key to access this api', type: 'string',
      },
    },
    required: ['email', 'oauth'],
    ...oauthLoginExample,
  },
  updateUser: {
    properties: {
      firstName: {
        description: 'user first name', type: 'string',
      },
      lastName: {
        description: 'user last name', type: 'string',
      },
      password: {
        description: 'password', type: 'string',
      },
      confirmPassword: {
        description: 'confirm password', type: 'string',
      },
    },
    required: ['firstName', 'lastName', 'password', 'confirmPassword'],
    ...updateUserExample,
  },
  addProject: {
    properties: {
      project: {
        name: {
          description: 'project name', type: 'string',
        },
        desc: {
          description: 'project description', type: 'string',
        },
      },
      pipeline: {
        description: 'pipeline type', type: 'string',
      },
    },
    required: ['project', 'pipeline'],
    ...addProjectExample,
  },
  updateProject: {
    properties: {
      name: {
        description: 'project name', type: 'string',
      },
      desc: {
        description: 'project description', type: 'string',
      },
      public: {
        description: 'is a public project?', type: 'boolean',
      },
      sharedTo: {
        description: 'other users', type: 'array',
      },
    },
    required: [],
    ...updateProjectExample,
  },
  projectAction: {
    properties: {
      code: {
        description: 'project unique id', type: 'string',
      },
    },
    required: ['code'],
    ...projectActionExample,
  },
  projectActionSuccessful: {
    properties: {
      project: {
        description: 'project json', type: 'object',
      },
      message: {
        description: 'action successful message', type: 'string',
      },
      success: {
        description: 'true or false', type: 'boolean',
      },
    },
    ...actionSuccessfulExample,
  },
  projectsActionSuccessful: {
    properties: {
      projects: {
        description: 'list of projects', type: 'array',
      },
      message: {
        description: 'action successful message', type: 'string',
      },
      success: {
        description: 'true or false', type: 'boolean',
      },
    },
    ...actionSuccessfulExample,
  },
  addUpload: {
    properties: {
      file: {
        description: 'File', type: 'string',
      },
      name: {
        description: 'File name', type: 'string',
      },
      desc: {
        description: 'File description', type: 'string',
      },
      type: {
        description: 'File type', type: 'string',
      },
      size: {
        description: 'File size', type: 'integer',
      },
    },
    required: ['name', 'type', 'size'],
    ...addUploadExample,
  },
  updateUpload: {
    properties: {
      name: {
        description: 'File name', type: 'string',
      },
      desc: {
        description: 'File description', type: 'string',
      },
      public: {
        description: 'is a public file?', type: 'boolean',
      },
      sharedTo: {
        description: 'other users', type: 'array',
      },
    },
    required: [],
    ...updateUploadExample,
  },
  projectFiles: {
    properties: {
      fileTypes: {
        description: 'File extensions', type: 'array',
      },
    },
    required: [],
    ...projectFilesExample,
  },
  actionSuccessful: {
    properties: {
      message: {
        description: 'action successful message', type: 'string',
      },
      success: {
        description: 'true or false', type: 'boolean',
      },
    },
    ...actionSuccessfulExample,
  },
  actionFailed: {
    properties: {
      message: {
        description: 'action failed message', type: 'string',
      },
      success: {
        description: 'true or false', type: 'boolean',
      },
    },
    ...actionFailedExample,
  },
  serverError: {
    properties: {
      message: {
        description: 'error message', type: 'string',
      },
      success: {
        description: 'true or false', type: 'boolean',
      },
    },
    ...serverErrorExample,
  },
};
