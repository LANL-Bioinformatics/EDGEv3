import React, { useState, useEffect } from 'react'
import { Col, Row, Input } from 'reactstrap'
import { useForm } from 'react-hook-form'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const RangeTextInput = (props) => {
  const componentName = 'rangeTextInput'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const {
    register,
    formState: { errors },
    watch,
    trigger,
  } = useForm({
    mode: 'onChange',
  })

  const startReg = register('start', {
    setValueAs: (v) => parseInt(v),
    min: { value: props.minStart, message: 'Value must be >= ' + props.minStart },
    max: { value: props.maxStart, message: 'Value must be <= ' + props.maxStart },
    validate: (value) => {
      if (!value) {
        return 'An integer is required'
      } else if (!/^[0-9]+$/.test(value)) {
        return 'Not an integer.'
      }
    },
  })

  const endReg = register('end', {
    setValueAs: (v) => parseInt(v),
    min: { value: props.minEnd, message: 'Value must be >= ' + props.minEnd },
    max: { value: props.maxEnd, message: 'Value must be <= ' + props.maxEnd },
    validate: (value) => {
      if (!value) {
        return 'An integer is required'
      } else if (!/^[0-9]+$/.test(value)) {
        return 'Not an integer.'
      }
      const max = value
      const min = watch('start')
      if (min) {
        return max >= min || 'Value is less than ' + props.startText
      }
    },
  })

  const setNewState = (e) => {
    setState({
      ...form,
      [e.target.name]: e.target.value,
    })
    setDoValidation(doValidation + 1)
  }

  const setNewState2 = (name, value) => {
    if (form[name] === value) {
      return
    }
    setState({
      ...form,
      [name]: value,
    })
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form.start = props.defaultValueStart
    form.end = props.defaultValueEnd
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  //trigger validation method when any input changes
  useEffect(() => {
    trigger().then((result) => {
      form.validForm = result
      if (result) {
        form.errMessage = ''
      } else {
        let errMessage = ''
        if (errors.start) {
          errMessage += errors.start.message + '\n'
        }
        if (errors.end) {
          errMessage += errors.end.message + '\n'
        }
        if (errMessage) {
          form.errMessage = errMessage
          form.validForm = false
        }
      }
      //force updating parent's inputParams
      props.setParams(form, props.name)
    })
  }, [doValidation])

  return (
    <>
      <Row>
        <Col xs="12" md="6">
          {props.startText}
          <Input
            type="number"
            name="start"
            defaultValue={props.defaultValueStart}
            placeholder={props.placeholderStart}
            style={errors['start'] ? defaults.inputStyleWarning : defaults.inputStyle}
            onChange={(e) => {
              startReg.onChange(e) // method from hook form register
              setNewState2(e.target.name, parseInt(e.target.value)) // your method
            }}
            innerRef={startReg.ref}
            onKeyPress={(e) => {
              e.key === 'Enter' && e.preventDefault()
            }}
          />
          {errors.start && <p className="edge-form-input-error">{errors.start.message}</p>}
        </Col>
        <Col xs="12" md="6">
          {props.endText}
          <Input
            type="number"
            name="end"
            defaultValue={props.defaultValueEnd}
            placeholder={props.placeholderEnd}
            style={errors['end'] ? defaults.inputStyleWarning : defaults.inputStyle}
            onChange={(e) => {
              endReg.onChange(e) // method from hook form register
              setNewState2(e.target.name, parseInt(e.target.value)) // your method
            }}
            innerRef={endReg.ref}
            onKeyPress={(e) => {
              e.key === 'Enter' && e.preventDefault()
            }}
          />
          {errors.end && <p className="edge-form-input-error">{errors.end.message}</p>}
        </Col>
      </Row>
      <br></br>
    </>
  )
}
