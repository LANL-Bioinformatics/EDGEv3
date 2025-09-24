import React, { useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Link, useNavigate, useLocation } from 'react-router-dom'
import { Button, Col, Row } from 'reactstrap'

import { login } from 'src/redux/reducers/edge/userSlice'
import { cleanMessage, cleanError } from '../../../redux/reducers/messageSlice'
import LoginForm from './forms/LoginForm'
import config from 'src/config'

const Login = (props) => {
  const navigate = useNavigate()
  const location = useLocation()
  const dispatch = useDispatch()
  const user = useSelector((state) => state.user)
  const loginErrors = useSelector((state) => state.message.errors)
  const messages = useSelector((state) => state.message.messages)

  useEffect(() => {
    dispatch(cleanMessage())
    dispatch(cleanError())
  }, [])

  useEffect(() => {
    if (user.isAuthenticated) {
      if (location.state) {
        navigate(location.state.from)
      } else {
        navigate('/home')
      }
    }
  }, [user, props, navigate, location])

  const handleValidSubmit = (data) => {
    const userData = { ...data }
    dispatch(login(userData))
  }

  const toOrcidLogin = () => {
    navigate('/oauth')
  }

  return (
    <div className="app-body flex-row align-items-center">
      <Row className="justify-content-center um-form">
        <Col md="4">
          <LoginForm messages={messages} errors={loginErrors} onSubmit={handleValidSubmit} />
          <Row className="justify-content-center">
            {config.APP.EMAIL_IS_ENABLED && (
              <>
                <Col xs="12">
                  <Link to="/activate">
                    <Button color="link" className="px-0">
                      Account not active?
                    </Button>
                  </Link>
                </Col>
                <Col xs="12">
                  <Link to="/resetPassword">
                    <Button color="link" className="px-0">
                      Forgot your password?
                    </Button>
                  </Link>
                </Col>
              </>
            )}
            <Col xs="12">
              <Link to="/register">
                <Button color="link" className="px-0">
                  No account? Create one
                </Button>
              </Link>
            </Col>
            {config.ORCID.IS_ENABLED && (
              <>
                <Col xs="12">
                  <hr></hr>
                  <p className="text-muted">
                    Use your ORCID account for faster login or registration
                  </p>

                  <Button
                    color="success"
                    className="px-4 text-white"
                    block
                    onClick={() => toOrcidLogin()}
                  >
                    Login with ORCiD
                  </Button>
                </Col>
              </>
            )}
          </Row>
        </Col>
      </Row>
    </div>
  )
}

export default Login
