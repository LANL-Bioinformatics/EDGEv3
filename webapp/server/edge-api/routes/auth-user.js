const router = require('express').Router();
const { validationRules: updateValidationRules, validate: updateValidate } = require('../validations/user-update-validator');
const {
  info,
  updateOne,
  getUsers,
} = require('../controllers/auth-user-controller');

/**
 * @swagger
 * /api/auth-user/info:
 *   get:
 *     summary: Find out if user can submit a new job
 *     tags: [AuthUser]
 *     security:
 *       - bearerAuth: []
 *     responses:
 *       200:
 *         description: Action successful.
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionSuccessful'
 *       400:
 *         description: Invalid input
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionFailed'
 *       500:
 *         description: API server error
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/serverError'
 */
router.get('/info', async (req, res) => {
  await info(req, res);
});

/**
 * @swagger
 * /api/auth-user:
 *   put:
 *     summary: Update current user
 *     tags: [AuthUser]
 *     security:
 *       - bearerAuth: []
 *     requestBody:
 *       required: true
 *       content:
 *         application/json:
 *           schema:
 *             $ref: '#/components/models/updateUser'
 *     responses:
 *       200:
 *         description: Action successful.
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionSuccessful'
 *       400:
 *         description: Invalid input
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionFailed'
 *       500:
 *         description: API server error
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/serverError'
 */
router.put('/', updateValidationRules(), updateValidate, async (req, res) => {
  await updateOne(req, res);
});

/**
 * @swagger
 * /api/auth-user/users:
 *   get:
 *     summary: Find all active users for share/unshare actions
 *     tags: [AuthUser]
 *     security:
 *       - bearerAuth: []
 *     responses:
 *       200:
 *         description: Action successful.
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionSuccessful'
 *       400:
 *         description: Invalid input
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/actionFailed'
 *       500:
 *         description: API server error
 *         content:
 *           application/json:
 *             schema:
 *               $ref: '#/components/models/serverError'
 */
router.get('/users', async (req, res) => {
  await getUsers(req, res);
});
module.exports = router;
