import React, { useState, useEffect, useMemo, Suspense } from 'react'
import { Col, Row, Button } from 'reactstrap'
import DropdownTreeSelect from 'react-dropdown-tree-select'
import 'react-dropdown-tree-select/dist/styles.css'
import { MyTooltip, ErrorTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const TreeSelectInput = (props) => {
  const componentName = 'treeSelectInput'
  const [resetTreeSelect, setResetTreeSelect] = useState(0)
  const [form, setState] = useState({ ...components[componentName] })

  const onSelectionChange = (currentNode, selectedNodes) => {
    if (!currentNode.label) {
      return
    }
    setState({ ...form, selections: selectedNodes })
  }

  const treeSelectorSearchPredicate = (node, searchTerm) => {
    //slow
    //only show child nodes: node._depth === 1
    return node.label && node.label.toLowerCase().includes(searchTerm) && node._depth === 1
    //fast
    //return node.label && node.label.toLowerCase().startsWith(searchTerm)
  }
  //Only re-rendered when data changes
  const treeSelector = useMemo(
    () => (
      <DropdownTreeSelect
        id={'treeSelect'}
        data={props.data}
        searchPredicate={treeSelectorSearchPredicate}
        className="edge-tree-select"
        texts={{ placeholder: props.placeholder }}
        clearSearchOnChange={false}
        keepOpenOnSelect={true}
        keepTreeOnSearch={true}
        keepChildrenOnSearch={false}
        showChildren={true}
        mode={props.mode}
        showPartiallySelected={false}
        inlineSearchInput={false}
        onChange={(currentNode, selectedNodes) => onSelectionChange(currentNode, selectedNodes)}
      />
    ),
    [props.data, props.reset, resetTreeSelect],
  )

  useEffect(() => {
    setState({ ...form, selections: props.value ? props.value : [] })
  }, [props.reset, resetTreeSelect]) // eslint-disable-line react-hooks/exhaustive-deps

  //trigger validation method when input changes
  useEffect(() => {
    form.validForm = true
    form.errMessage = ''
    const min = props.min ? props.min : 1
    const max = props.max ? props.max : 1000000
    if (!props.isOptional && (form.selections.length < min || form.selections.length > max)) {
      form.validForm = false
      form.errMessage = `Select at least ${min} but no more than ${max} items`
    }
    props.setParams(form, props.name)
  }, [form]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`textInputTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
            />
          ) : (
            <>
              {props.text}
              {errors && errors.textInput && props.showErrorTooltip && (
                <ErrorTooltip
                  id={`textInputErrTooltip-${props.name}`}
                  tooltip={errors.textInput.message}
                />
              )}
            </>
          )}
          <br></br>
          {props.showSelections && form.selections.length > 0 && (
            <>
              <span className="text-muted edge-text-size-small">
                {form.selections.length} {props.showSelectionsText}
              </span>
              <br></br>
              <Button
                size="sm"
                color="outline-primary"
                onClick={(e) => {
                  setResetTreeSelect(resetTreeSelect + 1)
                }}
              >
                Clear selection(s)
              </Button>
            </>
          )}
        </Col>
        <Col xs="12" md="9">
          {props.data && treeSelector}
        </Col>
      </Row>
    </>
  )
}
