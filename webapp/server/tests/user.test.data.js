module.exports = {
  userSuccess: {
    firstName: 'edge',
    lastName: 'lanl',
    email: 'my@edge.lanl',
    password: 'my#4EDGE',
    confirmPassword: 'my#4EDGE',
    active: false,
    // use 'test: true' to disable sendmail
    test: true,
  },
  registerInvalidInput: {
    firstName: 'not-good',
    lastName: '',
    email: 'my@edge',
    password: 'my#fEDGE',
    confirmPassword: 'my#5EDGE',
    active: 'bool',
  },
  activateInvalidInput: {
    email: 'notEmail',
    token: '',
  },
  activateUserNotFound: {
    email: 'notfound@edge.lanl',
    token: 'somestring',
  },
  getActionLinkInvalidInput: {
    email: 'notEmail',
    actionURL: 'not a URL',
  },
  getActionLinkNotFound: {
    email: 'notfound@edge.lanl',
    actionURL: 'http://my.edge/activate',
  },
  resetPasswordInvalidInput: {
    email: 'notEmail',
    token: '',
    newPassword: '',
  },
  resetPasswordUserNotFound: {
    email: 'notfound@edge.lanl',
    token: 'somestring',
    newPassword: 'my#4EDGE',
  },
  loginInvalidInput: {
    email: 'notEmail',
    password: '',
  },
  loginUserNotFound: {
    email: 'notfound@edge.lanl',
    password: 'my#4EDGE',
  },
  oauthLoginSuccess: {
    firstName: 'edge',
    lastName: 'lanl',
    email: '0000-1111-1111-1111@orcid.org',
    oauth: 'orcid',
  },
  oauthLoginInvalidInput: {
    email: 'notEmail',
  },
  oauthLoginUserNotActive: {
    firstName: 'edge',
    lastName: 'lanl',
    email: 'my@edge.lanl',
    oauth: 'orcid',
  },
};
