import React, { useState, useEffect } from 'react'
import { Button, Col, Row } from 'reactstrap'
import { useForm, useFieldArray, Controller } from 'react-hook-form'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'
import { RangeTextInput } from './RangeTextInput'

export const RangeTextInputArray = (props) => {
  const componentName = 'rangeTextInputArray'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)
  const { control } = useForm({
    mode: defaults['form_mode'],
  })

  const {
    fields: rangeTextInputFields,
    append: rangeTextInputAppend,
    remove: rangeTextInputRemove,
  } = useFieldArray({
    control,
    name: 'rangeTextInput',
  })

  const setRangeTextInputArray = (params, index) => {
    form.rangeInputs[index] = { ...params }
    setDoValidation(doValidation + 1)
  }

  //default 0 dataset
  useEffect(() => {
    setState({ ...form, rangeInputs: [] })
    setDoValidation(doValidation + 1)
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  //trigger validation method when input changes
  useEffect(() => {
    form.validForm = true
    form.errMessage = ''
    form.rangeInputs.forEach((item) => {
      if (!item.validForm) {
        form.validForm = false
        form.errMessage = 'Input error'
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
              id={`rangeTextInputArrayTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
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
              rangeTextInputAppend({ name: 'rangeTextInput' })
              setDoValidation(doValidation + 1)
            }}
          >
            Add more {props.rangeText}
          </Button>
        </Col>
      </Row>
      <br></br>
      {rangeTextInputFields.map((item, index) => (
        <div key={item.id}>
          <Row>
            <Col md="3" className="edge-sub-field">
              {props.rangeText} {index + 1}
            </Col>
            <Col xs="12" md="9">
              <Controller
                render={({ field: { ref, ...rest }, fieldState }) => (
                  <RangeTextInput
                    index={index}
                    name={index}
                    rangeText={props.rangeText}
                    setParams={setRangeTextInputArray}
                    startText={props.startText}
                    endText={props.endText}
                    defaultValueStart={props.defaultValueStart}
                    minStart={props.minStart}
                    maxStart={props.maxStart}
                    defaultValueEnd={props.defaultValueEnd}
                    minEnd={props.minEnd}
                    maxEnd={props.maxEnd}
                  />
                )}
                name={`rangeTextInput[${index}]`}
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
                  form.rangeInputs.splice(index, 1)
                  rangeTextInputRemove(index)
                  setDoValidation(doValidation + 1)
                }}
              >
                Remove {props.rangeText} {index + 1}
              </Button>
            </Col>
          </Row>
          <br></br>
        </div>
      ))}
    </>
  )
}
