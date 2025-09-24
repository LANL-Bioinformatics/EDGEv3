import React, { Suspense, useEffect } from 'react'
import { BrowserRouter as Router, Route, Routes } from 'react-router-dom'
import { useSelector } from 'react-redux'

import { CSpinner, useColorModes } from '@coreui/react'
import './scss/style.scss'

import { setUser, logout } from 'src/redux/reducers/edge/userSlice'
import { jwtDecode } from 'jwt-decode'
import { setAuthToken } from 'src/edge/common/util'
import store from 'src/redux/store'

// Check for token to keep user logged in
if (localStorage.jwtToken) {
  // Set auth token header auth
  const token = localStorage.jwtToken
  setAuthToken(token)
  // Decode token and get user info and exp
  const decoded = jwtDecode(token)
  // Set user and isAuthenticated
  store.dispatch(
    setUser({
      isAuthenticated: true,
      profile: decoded,
    }),
  )
} else {
  store.dispatch(logout())
}
//logout all tabs
window.addEventListener('storage', (e) => {
  if (e.key === 'jwtToken' && e.oldValue && !e.newValue) {
    store.dispatch(logout())
  }
})

// Containers
const DefaultLayout = React.lazy(() => import('./layout/DefaultLayout'))

const App = () => {
  const { isColorModeSet, setColorMode } = useColorModes('')
  const storedTheme = useSelector((state) => state.theme)

  useEffect(() => {
    setColorMode(storedTheme)
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <Router>
      <Suspense
        fallback={
          <div className="pt-3 text-center">
            <CSpinner color="primary" variant="grow" />
          </div>
        }
      >
        <Routes>
          <Route path="*" name="Home" element={<DefaultLayout />} />
        </Routes>
      </Suspense>
    </Router>
  )
}

export default App
