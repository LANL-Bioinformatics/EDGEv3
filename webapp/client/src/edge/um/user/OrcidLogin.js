import React, { useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { useNavigate, useLocation } from 'react-router-dom'
import { jwtDecode } from 'jwt-decode'
import queryString from 'query-string'
import { LoaderDialog, MessageDialog } from '../../common/Dialogs'
import { login } from 'src/redux/reducers/edge/userSlice'
import { cleanMessage } from 'src/redux/reducers/messageSlice'
import config from 'src/config'

const OrcidLogin = (props) => {
  const navigate = useNavigate()
  const location = useLocation()
  const dispatch = useDispatch()
  const errors = useSelector((state) => state.message.errors)
  const messages = useSelector((state) => state.message.messages)
  const page = useSelector((state) => state.page)
  const user = useSelector((state) => state.user)

  useEffect(() => {
    if (user.isAuthenticated) {
      console.log(localStorage.loginFrom)
      if (localStorage.loginFrom) {
        navigate(localStorage.loginFrom)
        // Remove from local storage
        localStorage.removeItem('loginFrom')
      } else {
        navigate('/home')
      }
    }
  }, [user])

  const getORCID = () => {
    const url = config.ORCID.AUTH_URI
    //console.log(url)
    //open ORCiD login page
    window.open(url, '_self')
  }
  const HandleMessages = (message, token) => {
    // Do stuff
    //console.log("got mes", message)
    const data = message
    const userData = {
      firstName: data.given_name,
      lastName: data.family_name,
      email: data.sub + '@orcid.org',
      oauth: 'orcid',
    }

    //clean up error messages
    dispatch(cleanMessage())
    //orcid login
    const status = 'active'
    if (!userData.lastName) {
      //is optional in ORCiD
      userData.lastName = 'unknown'
    }
    userData.status = status
    dispatch(login(userData))
  }

  const closeMsgModal = () => {
    //clean up error messages
    dispatch(cleanMessage())
  }
  useEffect(() => {
    if (location && location.hash) {
      const parsed = queryString.parse(location.hash)
      //console.log("parsed", parsed.id_token)
      if (parsed.id_token) {
        const decoded = jwtDecode(parsed.id_token)
        HandleMessages(decoded, parsed.id_token)
      }
    } else {
      getORCID()
    }
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <LoaderDialog loading={page.submitting_form} text="Verifying email..." />
      <MessageDialog
        isOpen={messages.sociallogin ? true : false}
        message={messages.sociallogin}
        handleClickClose={closeMsgModal}
      />
      <MessageDialog
        isOpen={errors.sociallogin ? true : false}
        message={errors.sociallogin}
        handleClickClose={closeMsgModal}
        className="modal-danger"
      />
    </>
  )
}

export default OrcidLogin
