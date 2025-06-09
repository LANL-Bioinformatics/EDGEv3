const router = require('express').Router();
const passport = require('passport');
const userRoutes = require('./edge-api/routes/user');
const projectRoutes = require('./edge-api/routes/project');
const authUserRoutes = require('./edge-api/routes/auth-user');
const authUserProjectRoutes = require('./edge-api/routes/auth-user-project');
const authUserUploadRoutes = require('./edge-api/routes/auth-user-upload');
const authUserDataRoutes = require('./edge-api/routes/auth-user-data');
const adminUserRoutes = require('./edge-api/routes/admin-user');
const adminProjectRoutes = require('./edge-api/routes/admin-project');
const adminUploadRoutes = require('./edge-api/routes/admin-upload');
const metagRoutes = require('./workflow_api/routes/metag');

/* GET home page. */
router.get('/', (req, res) => {
  res.send('The API server is on!');
});

router.use('/user', userRoutes);
router.use('/public/projects', projectRoutes);
router.use('/auth-user', passport.authenticate('user', { session: false }), authUserRoutes);
router.use('/auth-user', passport.authenticate('user', { session: false }), authUserProjectRoutes);
router.use('/auth-user', passport.authenticate('user', { session: false }), authUserUploadRoutes);
router.use('/auth-user', passport.authenticate('user', { session: false }), authUserDataRoutes);
router.use('/admin', passport.authenticate('admin', { session: false }), adminUserRoutes);
router.use('/admin', passport.authenticate('admin', { session: false }), adminProjectRoutes);
router.use('/admin', passport.authenticate('admin', { session: false }), adminUploadRoutes);
router.use('/workflow/metag', metagRoutes);

module.exports = router;
