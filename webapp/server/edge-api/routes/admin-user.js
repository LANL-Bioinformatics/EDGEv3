const router = require('express').Router();
const { validationRules: registerValidationRules, validate: registerValidate } = require('../validations/user-validator');
const { validationRules: updateValidationRules, validate: updateValidate } = require('../validations/user-update-validator');
const { validationRules: deleteValidationRules, validate: deleteValidate } = require('../validations/admin-user-validator');
const {
  addOne,
  updateOne,
  removeOne,
  getAll,
} = require('../controllers/admin-user-controller');

/**
 * @swagger
 * /api/admin/users:
 *   get:
 *     summary: Find all users
 *     tags: [Admin]
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
  await getAll(req, res);
});

/**
 * @swagger
 * /api/admin/users:
 *   post:
 *     summary: Create a new user
 *     tags: [Admin]
 *     security:
 *       - bearerAuth: []
 *     requestBody:
 *       required: true
 *       content:
 *         application/json:
 *           schema:
 *             $ref: '#/components/models/registerUser'
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
router.post('/users', registerValidationRules(), registerValidate, async (req, res) => {
  await addOne(req, res);
});

/**
 * @swagger
 * /api/admin/users/{user}:
 *   put:
 *     summary: Update existing user
 *     tags: [Admin]
 *     security:
 *       - bearerAuth: []
 *     parameters:
 *      - in: path
 *        name: user
 *        required: true
 *        type: string
 *        value: my@edge.lanl
 *        description: The user email
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
router.put('/users/:user', updateValidationRules(), updateValidate, async (req, res) => {
  await updateOne(req, res);
});

/**
 * @swagger
 * /api/admin/users/{user}:
 *   delete:
 *     summary: delete user if there is no project associated with the user
 *     tags: [Admin]
 *     security:
 *       - bearerAuth: []
 *     parameters:
 *      - in: path
 *        name: user
 *        required: true
 *        type: string
 *        value: my@edge.lanl
 *        description: The user email
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
router.delete('/users/:user', deleteValidationRules(), deleteValidate, async (req, res) => {
  await removeOne(req, res);
});

module.exports = router;
