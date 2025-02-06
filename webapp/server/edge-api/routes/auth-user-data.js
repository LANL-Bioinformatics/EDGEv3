const router = require('express').Router();
const {
  getPublicFiles,
  getGlobusFiles,
} = require('../controllers/auth-user-data-controller');

/**
 * @swagger
 * /api/auth-user/data/public:
 *   post:
 *     summary: List public data
 *     tags: [AuthUser]
 *     security:
 *       - bearerAuth: []
 *     requestBody:
 *       required: false
 *       content:
 *         application/json:
 *           schema:
 *             $ref: '#/components/models/publicData'
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

router.post('/data/public', async (req, res) => {
  await getPublicFiles(req, res);
});
/**
 * @swagger
 * /api/auth-user/data/globus:
 *   post:
 *     summary: List globus data
 *     tags: [AuthUser]
 *     security:
 *       - bearerAuth: []
 *     requestBody:
 *       required: false
 *       content:
 *         application/json:
 *           schema:
 *             $ref: '#/components/models/publicData'
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
router.post('/data/globus', async (req, res) => {
  await getGlobusFiles(req, res);
});

module.exports = router;
