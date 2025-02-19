import React, { useState, useEffect } from 'react'
import { Button, Col, Row } from 'reactstrap'
import { useForm, useFieldArray, Controller } from 'react-hook-form'

import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'
import { TextInput } from './TextInput'

export const TextInputArray = (props) => {
  const componentName = 'textInputArray'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)
  const { control } = useForm({
    mode: defaults['form_mode'],
  })

  const {
    fields: inputFields,
    append: inputAppend,
    remove: inputRemove,
  } = useFieldArray({
    control,
    name: 'input',
  })

  const setInputArray = (params, index) => {
    form.textInputs[index] = { ...params }
    setDoValidation(doValidation + 1)
  }

  //default 0 dataset
  useEffect(() => {
    setState({ ...form, ['textInputs']: [] })
    setDoValidation(doValidation + 1)
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  //trigger validation method when input changes
  useEffect(() => {
    form.validForm = true
    form.errMessage = ''
    form.textInputs.forEach((item) => {
      if (!item.validForm) {
        form.validForm = false
        form.errMessage = item.errMessage
      }
    })
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`textInputArrayTooltip-${props.name}`}
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
          <Button
            size="sm"
            className="btn-pill"
            color="outline-primary"
            onClick={() => {
              inputAppend({ name: 'input' })
              setDoValidation(doValidation + 1)
            }}
          >
            Add more {props.inputText}
          </Button>
        </Col>
      </Row>
      <br></br>
      {inputFields.map((item, index) => (
        <div key={item.id}>
          <Row>
            <Col md="3" className="edge-sub-field">
              {props.inputText} {index + 1}
            </Col>
            <Col xs="12" md="9">
              <Controller
                render={({ field: { ref, ...rest }, fieldState }) => (
                  <TextInput
                    index={index}
                    name={index}
                    defaultValue={''}
                    setParams={setInputArray}
                    isOptional={false}
                    isValidTextInput={props.isValidTextInput}
                    errMessage={'Required'}
                  />
                )}
                name={`input[${index}]`}
                control={control}
              />
            </Col>
          </Row>
          <Row>
            <Col md="3"></Col>
            <Col xs="12" md="9">
              <Button
                size="sm"
                className="btn-pill"
                color="ghost-primary"
                onClick={() => {
                  form.textInputs.splice(index, 1)
                  inputRemove(index)
                  setDoValidation(doValidation + 1)
                }}
              >
                Remove {props.inputText} {index + 1}
              </Button>
            </Col>
          </Row>
          <br></br>
        </div>
      ))}
    </>
  )
}
