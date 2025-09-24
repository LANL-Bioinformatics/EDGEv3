import React, { useState, useEffect } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Col, Row } from 'reactstrap'
import { ToastContainer } from 'react-toastify'
import 'react-toastify/dist/ReactToastify.css'

import { LoaderDialog } from '../../common/Dialogs'
import { notify } from '../../common/util'
import { update } from 'src/redux/reducers/edge/userSlice'
import ProfileForm from './forms/ProfileForm'

const Profile = (props) => {
  const dispatch = useDispatch()
  const user = useSelector((state) => state.user)
  const updateErrors = useSelector((state) => state.message.errors)
  const page = useSelector((state) => state.page)
  const [sendNotify, setSendNotify] = useState(false)

  const handleValidSubmit = async (data) => {
    setSendNotify(false)
    const newUser = {
      email: data.email,
      firstName: data.firstName,
      lastName: data.lastName,
      notification: data.notification,
    }
    if (data.notification) {
      newUser['notificationEmail'] = data.notificationEmail
    }
    if (data.changePw) {
      newUser['password'] = data.password
      newUser['confirmPassord'] = data.confirmPassord
    }
    dispatch(update(newUser))

    //wait for updating complete
    setTimeout(() => {
      setSendNotify(true)
    }, 200)
  }

  const notifyUpdateResult = () => {
    if (updateErrors['update']) {
      notify('error', 'Failed to update profile.')
    } else {
      notify('success', 'Your profile has been updated successfully.')
    }
  }

  useEffect(() => {
    if (sendNotify) {
      notifyUpdateResult()
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sendNotify])

  return (
    <div className="app-body flex-row align-items-center">
      <LoaderDialog loading={page.submitting_form} />
      <ToastContainer />
      <Row className="justify-content-center um-form">
        <Col md="4">
          {updateErrors.update && <p className="edge-form-input-error">{updateErrors.update}</p>}
          <ProfileForm user={user} loading={page.submitting_form} onSubmit={handleValidSubmit} />
        </Col>
      </Row>
    </div>
  )
}

export default Profile
