import React from 'react'
import { Form, Input, Button, Col, InputGroup, InputGroupText, Row } from 'reactstrap'

import { useForm } from 'react-hook-form'

import CIcon from '@coreui/icons-react'
import { cilUser, cilLockLocked } from '@coreui/icons'
import { MyTooltip } from '../../../common/MyTooltip'
import { umTips } from '../../defaults'

const LoginForm = (props) => {
  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm({
    mode: 'onChange',
  })

  const passwordReg = {
    ...register('password', {
      required: ' Password is required',
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

  const emailReg = {
    ...register('email', {
      required: 'Email is required',
      pattern: {
        // Validation pattern
        value: /^[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}$/i,
        message: 'Invalid email address',
      },
    }),
  }

  return (
    <Form onSubmit={handleSubmit(props.onSubmit)}>
      <span className="red-text">{props.messages.login}</span>
      <span className="red-text">{props.errors.login}</span>
      <h1>Login</h1>
      <p className="text-muted">Sign In to your account</p>
      <span className="red-text">{props.errors.email}</span>
      <InputGroup className="mb-3">
        <InputGroupText>
          <CIcon icon={cilUser} />
        </InputGroupText>
        <Input
          type="text"
          name="email"
          placeholder="Email"
          onChange={(e) => {
            emailReg.onChange(e) // method from hook form register
          }}
          innerRef={emailReg.ref}
        />
      </InputGroup>
      {errors.email && <p className="edge-form-input-error">{errors.email.message}</p>}
      <span className="red-text">{props.errors.password}</span>
      <InputGroup className="mb-4">
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
      <Row>
        <Col xs="12">
          <Button color="primary" type="submit" className="px-4" block>
            Login
          </Button>
        </Col>
      </Row>
    </Form>
  )
}

export default LoginForm
