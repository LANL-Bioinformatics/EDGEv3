const request = require('supertest');
const db = require('./db');
const app = require('./testServer');
const testData = require('./user.test.data');

// Setup connection to the database
beforeAll(async () => db.connect());
beforeEach(async () => db.clear());
afterAll(async () => db.close());

// Pass supertest agent for each test
const agent = request.agent(app);

// user test suite
describe('Test user APIs', () => {
  describe('Register', () => {
    it('/api/user/register successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      expect(response.status).toEqual(200);
      expect(response.body.user.notification.email).toEqual(testData.userSuccess.email);
    });
    it('/api/user/register failed: Email already exists', async () => {
      await agent.post('/api/user/register').send(testData.userSuccess);
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Email my@edge.lanl already exists');
    });
    it('/api/user/register failed: Invaid input', async () => {
      const response = await agent.post('/api/user/register').send(testData.registerInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.firstName).toEqual('Invalid firstName.');
      expect(response.body.error.lastName).toEqual('Invalid lastName.');
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.password).toEqual('Password must be at least 8 characters long and contain at least one uppercase, at least one lower case and at least one special character.');
      expect(response.body.error.confirmPassword).toEqual('Password confirmation does not match with password.');
      expect(response.body.error.active).toEqual('Invalid active. Must be true or false.');
    });
  });
  describe('Activate', () => {
    it('/api/user/activate successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      const response2 = await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      expect(response2.status).toEqual(200);
    });
    it('/api/user/activate failed: Invalid input', async () => {
      const response = await agent.put('/api/user/activate').send(testData.validInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.token).toEqual('Token should not be empty.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/activate failed: Account not found', async () => {
      const response = await agent.put('/api/user/activate').send(testData.activateUserNotFound);
      expect(response.status).toEqual(400);
      expect(response.body.error.activate).toEqual('Account notfound@edge.lanl not found');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/activate failed: Invalid token', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      const response2 = await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: 'wrongpassword',
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.activate).toEqual('Invalid token.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
    it('/api/user/activate failed: Account already activated', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      const response2 = await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.activate).toEqual('Your account has already been activated.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
  });
  describe('getActivationLink', () => {
    it('/api/user/getActivationLink successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      // use 'test: true' to disable sendmail
      const response2 = await agent.post('/api/user/getActivationLink').send({
        email: response.body.user.email, actionURL: 'http://my.edge/activate', test: true,
      });
      expect(response2.status).toEqual(200);
    });
    it('/api/user/getActivationLink failed: Invalid input', async () => {
      const response = await agent.post('/api/user/getActivationLink').send(testData.getActionLinkInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.actionURL).toEqual('Invalid action URL.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/getActivationLink failed: Account not found', async () => {
      const response = await agent.post('/api/user/getActivationLink').send(testData.getActionLinkNotFound);
      expect(response.status).toEqual(400);
      expect(response.body.error.getActivationLink).toEqual('Account notfound@edge.lanl not found');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/getActivationLink failed: Account already activated ', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      const response2 = await agent.post('/api/user/getActivationLink').send({
        email: response.body.user.email, actionURL: 'http://my.edge/activate',
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.getActivationLink).toEqual('Your account has already been activated.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
  });
  describe('ResetPassword', () => {
    it('/api/user/resetPassword successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      const response2 = await agent.put('/api/user/resetPassword').send({
        email: response.body.user.email, token: response.body.user.password, newPassword: 'new#4EDGE',
      });
      expect(response2.status).toEqual(200);
    });
    it('/api/user/resetPassword failed: Invalid input', async () => {
      const response = await agent.put('/api/user/resetPassword').send(testData.resetPasswordInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.token).toEqual('Token should not be empty.');
      expect(response.body.error.newPassword).toEqual('Password must be at least 8 characters long and contain at least one uppercase, at least one lower case and at least one special character.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/resetPassword failed: Account not found', async () => {
      const response = await agent.put('/api/user/resetPassword').send(testData.resetPasswordUserNotFound);
      expect(response.status).toEqual(400);
      expect(response.body.error.resetPassword).toEqual('Account notfound@edge.lanl not found');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/resetPassword failed: Account is not activated', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      const response2 = await agent.put('/api/user/resetPassword').send({
        email: response.body.user.email, token: response.body.user.password, newPassword: 'new#4EDGE',
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.resetPassword).toEqual('Your account is not active.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
    it('/api/user/resetPassword failed: Invalid token', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      const response2 = await agent.put('/api/user/resetPassword').send({
        email: response.body.user.email, token: 'wrongpassword', newPassword: 'new#4EDGE',
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.resetPassword).toEqual('Invalid token.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
  });
  describe('getResetPasswordLink', () => {
    it('/api/user/getResetPasswordLink successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      // use 'test: true' to disable sendmail
      const response2 = await agent.post('/api/user/getResetPasswordLink').send({
        email: response.body.user.email, actionURL: 'http://my.edge/activate', test: true,
      });
      expect(response2.status).toEqual(200);
    });
    it('/api/user/getResetPasswordLink failed: Invalid input', async () => {
      const response = await agent.post('/api/user/getResetPasswordLink').send(testData.getActionLinkInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.actionURL).toEqual('Invalid action URL.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/getResetPasswordLink failed: Account not found', async () => {
      const response = await agent.post('/api/user/getResetPasswordLink').send(testData.getActionLinkNotFound);
      expect(response.status).toEqual(400);
      expect(response.body.error.getResetPasswordLink).toEqual('Account notfound@edge.lanl not found');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/getResetPasswordLink failed: Account is not active ', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      const response2 = await agent.post('/api/user/getResetPasswordLink').send({
        email: response.body.user.email, actionURL: 'http://my.edge/activate',
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.getResetPasswordLink).toEqual('Your account is not active.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
  });
  describe('login', () => {
    it('/api/user/login successful', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      await agent.put('/api/user/activate').send({
        email: response.body.user.email, token: response.body.user.password,
      });
      const response2 = await agent.post('/api/user/login').send({
        email: response.body.user.email, password: testData.userSuccess.password,
      });
      expect(response2.status).toEqual(200);
    });
    it('/api/user/login failed: Invalid input', async () => {
      const response = await agent.post('/api/user/login').send(testData.loginInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.password).toEqual('Password must be at least 8 characters long and contain at least one uppercase, at least one lower case and at least one special character.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/login failed: Account not found', async () => {
      const response = await agent.post('/api/user/login').send(testData.loginUserNotFound);
      expect(response.status).toEqual(400);
      expect(response.body.error.login).toEqual('Account notfound@edge.lanl not found');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/login failed: Account is not active ', async () => {
      const response = await agent.post('/api/user/register').send(testData.userSuccess);
      const response2 = await agent.post('/api/user/login').send({
        email: response.body.user.email, password: testData.userSuccess.password,
      });
      expect(response2.status).toEqual(400);
      expect(response2.body.error.login).toEqual('Your account is not active.');
      expect(response2.body.message).toEqual('Action failed');
      expect(response2.body.success).toEqual(false);
    });
  });
  describe('oauthLogin', () => {
    it('/api/user/oauthLogin successful', async () => {
      const response = await agent.post('/api/user/oauthLogin').send(testData.oauthLoginSuccess);
      expect(response.status).toEqual(200);
    });
    it('/api/user/oauthLogin failed: Invalid input', async () => {
      const response = await agent.post('/api/user/oauthLogin').send(testData.oauthLoginInvalidInput);
      expect(response.status).toEqual(400);
      expect(response.body.error.firstName).toEqual('Invalid firstName.');
      expect(response.body.error.lastName).toEqual('Invalid lastName.');
      expect(response.body.error.email).toEqual('Invalid email address.');
      expect(response.body.error.oauth).toEqual('Oauth should not be empty.');
      expect(response.body.message).toEqual('Validation failed');
      expect(response.body.success).toEqual(false);
    });
    it('/api/user/oauthLogin failed: Account is not active ', async () => {
      await agent.post('/api/user/register').send(testData.userSuccess);
      const response = await agent.post('/api/user/oauthLogin').send(testData.oauthLoginUserNotActive);
      expect(response.status).toEqual(400);
      expect(response.body.error.oauthLogin).toEqual('Your account is not active.');
      expect(response.body.message).toEqual('Action failed');
      expect(response.body.success).toEqual(false);
    });
  });
});
