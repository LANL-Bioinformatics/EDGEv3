const fs = require('fs');
const Project = require('../models/project');
const Job = require('../models/job');
const { getAllFiles } = require('../../utils/common');
const { generateRunStats } = require('../../utils/cromwell');
const { generateRunStats: generateNextflowRunStats } = require('../../utils/nextflow');
const { generateWorkflowResult } = require('../../utils/workflow');
const config = require('../../config');

const getProject = async (code, type, user) => {
  try {
    // Use $eq to prevent query selector injections
    const project = await Project.findOne({ 'status': { $ne: 'delete' }, 'code': { $eq: code } });
    if (project === null) {
      return null;
    }
    if (type === 'admin') {
      return project;
    }
    if (type === 'user' && (project.owner === user.email || project.sharedTo.includes(user.email) || project.public)) {
      return project;
    } if (type === 'public' && project.public) {
      return project;
    } if (project.parent) {
      let query = { status: { $ne: 'delete' }, code: project.parent, public: true };
      if (type === 'user') {
        query = { status: { $ne: 'delete' }, code: project.parent, $or: [{ owner: user.email }, { sharedTo: user.email }, { public: true }] };
      }
      const pp = await Project.findOne(query);
      if (pp) {
        return project;
      }
      return null;
    }
    return null;
  } catch (err) {
    return Promise.reject(err);
  }
};

const updateProject = async (query, req) => {
  try {
    const proj = await Project.findOne(query);
    if (!proj) {
      return null;
    }
    if (req.body.name) {
      proj.name = req.body.name;
    }
    if (req.body.desc) {
      proj.desc = req.body.desc;
    }
    if (req.body.status) {
      proj.status = req.body.status;
    }
    if ('public' in req.body) {
      proj.public = req.body.public;
    }
    if (req.body.sharedTo) {
      proj.sharedTo = req.body.sharedTo;
    }
    if (req.body.jobPriority) {
      proj.jobPriority = req.body.jobPriority;
    }

    proj.updated = Date.now();
    const project = await proj.save();
    // set children to 'delete'
    if (project.type.match(/batch/) && project.status === 'delete') {
      await Project.find({
        status: { $ne: 'delete' }, code: { $in: project.children },
      }).then((projects) => {
        projects.forEach((subproject) => {
          const subproj = subproject;
          subproj.status = 'delete';
          subproj.updated = Date.now();
          subproj.save();
        });
      });
    }
    return project;
  } catch (err) {
    return Promise.reject(err);
  }
};

// Return conf.json as JSON object
const getProjectConf = async (code, type, req) => {
  try {
    const proj = await getProject(code, type, req.user);
    let confJson = {};
    if (proj) {
      const projHome = `${config.IO.PROJECT_BASE_DIR}/${code}`;
      confJson = JSON.parse(fs.readFileSync(`${projHome}/conf.json`));
    }
    return confJson;
  } catch (err) {
    return Promise.reject(err);
  }
};

const getProjectOutputs = async (code, type, req) => {
  try {
    const projDir = config.IO.PROJECT_BASE_DIR;
    const proj = await getProject(code, type, req.user);
    let files = [];
    if (proj && fs.existsSync(`${projDir}/${proj.code}/output`)) {
      files = getAllFiles(`${projDir}/${proj.code}/output`, files, req.body.fileTypes, '', `/projects/${proj.code}/output`, `${projDir}/${proj.code}/output`);
    }
    return files;
  } catch (err) {
    return Promise.reject(err);
  }
};

const getProjectBatchOutputs = (code, type, req) => new Promise((resolve, reject) => {
  const projDir = config.IO.PROJECT_BASE_DIR;
  getProject(code, type, req.user).then(project => {
    if (!project.type.includes('batch')) {
      reject(new Error('Not a batch project'));
    }
    Project.find({
      'status': { $ne: 'delete' }, 'code': { $in: project.children }
    }).then((projects) => {
      const allfiles = [];
      if (projects.length === 0) {
        resolve(allfiles);
      } else {
        projects.forEach((proj, index) => {
          const files = [];
          if (proj && fs.existsSync(`${projDir}/${proj.code}/output`)) {
            allfiles.push(...getAllFiles(`${projDir}/${proj.code}/output`, files, req.body.fileTypes, proj.name, `/projects/${proj.code}/output`, `${projDir}/${proj.code}/output`));
          }
          if (index === projects.length - 1) {
            resolve(allfiles);
          }
        });
      }
    });
  }).catch(err => { reject(err); });
});

const getProjectResult = async (code, type, req) => {
  try {
    const proj = await getProject(code, type, req.user);
    let result = {};
    if (proj) {
      const projHome = `${config.IO.PROJECT_BASE_DIR}/${code}`;
      const resultJson = `${projHome}/result.json`;

      if (!fs.existsSync(resultJson)) {
        generateWorkflowResult(proj);
      }
      result = JSON.parse(fs.readFileSync(resultJson));
    }
    return result;
  } catch (err) {
    return Promise.reject(err);
  }
};

const getProjectRunStats = async (code, type, req) => {
  try {
    const proj = await getProject(code, type, req.user);
    let stats = {};
    if (proj) {
      // get associated job
      const job = await Job.findOne({ 'project': { $eq: code } });
      if (!job) {
        return {};
      }
      const projHome = `${config.IO.PROJECT_BASE_DIR}/${code}`;
      const statsJson = `${projHome}/run_stats.json`;
      if (job.queue === 'nextflow') {
        await generateNextflowRunStats(proj);
      } else {
        generateRunStats(proj);
      }
      stats = JSON.parse(fs.readFileSync(statsJson));

    }
    return stats;
  } catch (err) {
    return Promise.reject(err);
  }
};

module.exports = {
  getAllFiles,
  getProject,
  updateProject,
  getProjectConf,
  getProjectOutputs,
  getProjectBatchOutputs,
  getProjectResult,
  getProjectRunStats,
};
