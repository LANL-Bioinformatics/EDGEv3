import React from 'react'

const Home = React.lazy(() => import('src/edge/Home'))
const PublicProjects = React.lazy(() => import('src/edge/um/public/Projects'))
const PublicProjectPage = React.lazy(() => import('src/workflows/project/projectPage/Public'))
const UserLogin = React.lazy(() => import('src/edge/um/user/Login'))
const OAuth = React.lazy(() => import('src/edge/um/user/OrcidLogin'))
// const UserRegister = React.lazy(() => import('src/edge/um/user/Register'))
// const UserActivate = React.lazy(() => import('src/edge/um/user/Activate'))
// const UserResetPassword = React.lazy(() => import('src/edge/um/user/ResetPassword'))

const routes = [
  { path: '/', exact: true, name: 'Home' },
  { path: '/home', name: 'Home', element: Home },
  { path: '/public/projects', name: 'PublicProjects', element: PublicProjects },
  { path: '/public/project', name: 'PublicProjectPage', element: PublicProjectPage },
  // admin login
  { path: '/edgeadminlogin', exact: true, name: 'Login', element: UserLogin },
  { path: '/oauth', name: 'OAuth', element: OAuth },
  // user/password login
  // { path: '/register', exact: true, name: 'Register', element: UserRegister },
  // { path: '/login', exact: true, name: 'Login', element: UserLogin },
  // { path: '/activate', exact: true, name: 'Activate', element: UserActivate },
  // { path: '/resetPassword', exact: true, name: 'ResetPassword', element: UserResetPassword },
]

export default routes
