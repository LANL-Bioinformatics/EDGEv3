const router = require('express').Router();
const passport = require('passport');
const userRoutes = require('./user');
const projectRoutes = require('./project');
const authUserRoutes = require('./auth-user');
const authUserProjectRoutes = require('./auth-user-project');
const authUserUploadRoutes = require('./auth-user-upload');
const authUserDataRoutes = require('./auth-user-data');
const adminUserRoutes = require('./admin-user');
const adminProjectRoutes = require('./admin-project');
const adminUploadRoutes = require('./admin-upload');

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

module.exports = router;
