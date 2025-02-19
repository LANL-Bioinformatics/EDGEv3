import React, { useEffect, useState } from 'react'
import { Form, Button, ButtonGroup, InputGroup, InputGroupText, Input } from 'reactstrap'
import { useForm } from 'react-hook-form'
import CIcon from '@coreui/icons-react'
import { cilUser, cilLockLocked } from '@coreui/icons'
import { MyTooltip } from '../../../common/MyTooltip'
import { umTips } from '../../defaults'
import config from 'src/config'

const Profile = (props) => {
  const [changePw, setChangePw] = useState(false)
  const [notification, setNotification] = useState(props.user.profile.notification.isOn)
  const isORCiDUser = props.user.profile.email.match(/\d+-\d+-\d+-\d+@orcid\.org/) ? true : false

  const {
    register,
    handleSubmit,
    formState: { errors },
    watch,
    setValue,
  } = useForm({
    mode: 'onChange',
  })

  const firstNameReg = {
    ...register('firstName', {
      required: 'Please enter your first name',
      pattern: { value: /^[A-Za-z]+$/, message: 'Your name must be composed only with letters' },
    }),
  }
  const lastNameReg = {
    ...register('lastName', {
      required: 'Please enter your last name',
      pattern: { value: /^[A-Za-z]+$/, message: 'Your name must be composed only with letters' },
    }),
  }

  const emailReg = {
    ...register('notificationEmail', {
      required: 'Email is required',
      pattern: {
        // Validation pattern
        value: /^[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}$/i,
        message: 'Invalid email address',
      },
    }),
  }

  const passwordReg = {
    ...register('password', {
      required: 'Please enter a password',
      minLength: { value: 8, message: 'Must be at least 8 characters long' },
      validate: {
        hasUpperCase: (value) =>
          /[A-Z]/.test(value) || 'Must contain at lease one uppercase letter',
        hasLowerCase: (value) =>
          /[a-z]/.test(value) || 'Must contain at lease one lowercase letter',
        hasNumber: (value) => /[0-9]/.test(value) || 'Must contain at lease one number',
        hasSpecialChar: (value) =>
          /[^A-Za-z0-9 ]/.test(value) || 'Must contain at lease one special character',
      },
    }),
  }

  const password2Reg = {
    ...register('confirmPassword', {
      validate: (value) => value === watch('password') || 'The passwords do not match',
    }),
  }

  useEffect(() => {
    setValue('firstName', props.user.profile.firstName)
    setValue('lastName', props.user.profile.lastName)
    setValue('email', props.user.profile.email)
    setValue(
      'notificationEmail',
      props.user.profile.notification.email
        ? props.user.profile.notification.email
        : props.user.profile.email,
    )
    setValue('notification', props.user.profile.notification.isOn)
    setNotification(props.user.profile.notification.isOn)
  }, [props.user, setValue])

  useEffect(() => {
    setValue('notification', notification)
    // disable form errors
    if (!notification) {
      setValue('notificationEmail', 'email@your.com')
    } else {
      setValue(
        'notificationEmail',
        props.user.profile.notification.email
          ? props.user.profile.notification.email
          : props.user.profile.email,
      )
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [notification])

  useEffect(() => {
    setValue('changePw', changePw)
    // disable form errors
    if (!changePw) {
      setValue('password', 'somePW#2')
      setValue('confirmPassword', 'somePW#2')
    } else {
      setValue('password', '')
      setValue('confirmPassword', '')
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [changePw])

  return (
    <Form onSubmit={handleSubmit(props.onSubmit)}>
      <h4 className="pt-3">Profile</h4>
      <hr></hr>

      <InputGroup className="mb-3">
        <InputGroupText>
          <CIcon icon={cilUser} />
        </InputGroupText>
        <Input type="text" name="email" placeholder={props.user.profile.email} disabled />
      </InputGroup>

      <InputGroup className="mb-3">
        <InputGroupText> First Name </InputGroupText>
        <Input
          type="text"
          name="firstName"
          onChange={(e) => {
            firstNameReg.onChange(e) // method from hook form register
          }}
          innerRef={firstNameReg.ref}
        />
      </InputGroup>
      {errors.firstName && <p className="edge-form-input-error">{errors.firstName.message}</p>}
      <InputGroup className="mb-3">
        <InputGroupText> Last Name </InputGroupText>
        <Input
          type="text"
          name="lastName"
          onChange={(e) => {
            lastNameReg.onChange(e) // method from hook form register
          }}
          innerRef={lastNameReg.ref}
        />
      </InputGroup>
      {errors.lastName && <p className="edge-form-input-error">{errors.lastName.message}</p>}
      <br></br>
      {config.APP.EMAIL_IS_ENABLED && (
        <>
          <b>Project Status Notification</b>
          <br></br>
          <div className="d-sm-inline-block">
            <ButtonGroup className="mr-3" aria-label="First group" size="sm">
              <Button
                color="outline-primary"
                onClick={() => {
                  setNotification(true)
                }}
                active={notification}
              >
                On
              </Button>
              <Button
                color="outline-primary"
                onClick={() => {
                  setNotification(false)
                }}
                active={!notification}
              >
                Off
              </Button>
            </ButtonGroup>
          </div>
          <br />
          <br />
          <Input type="hidden" name="notification" value={notification} />
          {notification && (
            <>
              <InputGroup className="mb-3">
                <InputGroupText>Email</InputGroupText>
                <Input
                  type="text"
                  name="notificationEmail"
                  onChange={(e) => {
                    emailReg.onChange(e) // method from hook form register
                  }}
                  innerRef={emailReg.ref}
                />
              </InputGroup>
              {errors.notificationEmail && (
                <p className="edge-form-input-error">{errors.notificationEmail.message}</p>
              )}
            </>
          )}
          <br></br>
        </>
      )}
      {!isORCiDUser && (
        <>
          <b>Change Password</b>
          <br></br>
          <div className="d-sm-inline-block">
            <ButtonGroup className="mr-3" aria-label="First group" size="sm">
              <Button
                color="outline-primary"
                onClick={() => {
                  setChangePw(true)
                }}
                active={changePw}
              >
                Yes
              </Button>
              <Button
                color="outline-primary"
                onClick={() => {
                  setChangePw(false)
                }}
                active={!changePw}
              >
                No
              </Button>
            </ButtonGroup>
          </div>
          <br />
          <br />
        </>
      )}
      <Input type="hidden" name="changePw" value={changePw} />
      {changePw && (
        <div>
          <InputGroup className="mb-3">
            <InputGroupText>
              <CIcon icon={cilLockLocked} />
            </InputGroupText>
            <Input
              type="password"
              name="password"
              placeholder="Password"
              onChange={(e) => {
                passwordReg.onChange(e) // method from hook form register
              }}
              innerRef={passwordReg.ref}
            />
          </InputGroup>
          {errors.password && (
            <p className="edge-form-input-error">
              {errors.password.message}
              <MyTooltip id="passwordRegister" tooltip={umTips.passwordHints} showTooltip={true} />
            </p>
          )}
          <InputGroup className="mb-4">
            <InputGroupText>
              <CIcon icon={cilLockLocked} />
            </InputGroupText>
            <Input
              type="password"
              name="confirmPassword"
              placeholder="Repeat password"
              onChange={(e) => {
                password2Reg.onChange(e) // method from hook form register
              }}
              innerRef={password2Reg.ref}
            />
          </InputGroup>
          {errors.confirmPassword && (
            <p className="edge-form-input-error">{errors.confirmPassword.message}</p>
          )}
        </div>
      )}
      <Button color="primary" disabled={props.loading} type="submit" block>
        Save Changes
      </Button>
    </Form>
  )
}

export default Profile
