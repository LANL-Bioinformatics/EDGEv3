import React, { useState, useEffect } from 'react'
import { Col, Row, Input } from 'reactstrap'
import { useForm } from 'react-hook-form'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const IntegerInput = (props) => {
  const componentName = 'integerInput'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const {
    register,
    formState: { errors },
    trigger,
  } = useForm({
    mode: defaults['form_mode'],
  })

  const integerInputReg = register('integerInput', {
    required: 'required, an integer',
    setValueAs: (v) => parseInt(v),
    min: { value: props.min, message: 'Value is less than ' + props.min },
    max: { value: props.max, message: 'Value is greater than ' + props.max },
    validate: (value) => {
      if (!/^[0-9]+$/.test(value)) {
        return 'Not an integer.'
      }
    },
  })

  const setNewState2 = (name, value) => {
    setState({
      ...form,
      [name]: value,
    })
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form.integerInput = props.defaultValue
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    //validate form
    trigger().then((result) => {
      form.validForm = result
      if (result) {
        form.errMessage = ''
      } else {
        let errMessage = props.name + ' input error'
        form.errMessage = errMessage
      }
      //force updating parent's inputParams
      props.setParams(form, props.name)
    })
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`integerInputTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
            />
          ) : (
            <>{props.text}</>
          )}
        </Col>
        <Col xs="12" md="9">
          <Input
            type="number"
            name="integerInput"
            key={props.name}
            defaultValue={props.defaultValue}
            style={errors['integerInput'] ? defaults.inputStyleWarning : defaults.inputStyle}
            onChange={(e) => {
              integerInputReg.onChange(e) // method from hook form register
              setNewState2(e.target.name, parseInt(e.target.value)) // your method
            }}
            innerRef={integerInputReg.ref}
          />
        </Col>
      </Row>
    </>
  )
}
