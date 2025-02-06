import React, { useState, useEffect } from 'react'
import { Button, Form } from 'reactstrap'
import { ToastContainer } from 'react-toastify'
import 'react-toastify/dist/ReactToastify.css'

import { workflowList } from 'src/util'
import { postData, getData, notify, apis } from '../../common/util'
import { LoaderDialog, MessageDialog } from '../../common/Dialogs'
import { HtmlText } from '../../common/HtmlText'
import SraDataTable from '../../um/common/SraDataTable'
import { Sra2fastq } from './forms/Sra2fastq'
import { workflowOptions } from './defaults'

const Main = (props) => {
  const [submitting, setSubmitting] = useState(false)
  const [requestSubmit, setRequestSubmit] = useState(false)
  const [selectedWorkflows, setSelectedWorkflows] = useState({})
  const [doValidation, setDoValidation] = useState(0)
  const [workflow, setWorkflow] = useState(workflowOptions[0].value)
  const [resetText, setResetText] = useState(0)
  const [openDialog, setOpenDialog] = useState(false)
  const [disabled, setDisabled] = useState(false)
  const [sysMsg, setSysMsg] = useState()

  const setWorkflowParams = (params, workflowName) => {
    //console.log('workflow:', params, workflowName)
    setSelectedWorkflows({ ...selectedWorkflows, [workflowName]: params })
    setDoValidation(doValidation + 1)
  }

  //submit button clicked
  const onSubmit = () => {
    setSubmitting(true)
    let formData = {}
    formData.category = workflowList[workflow].category
    // set project info
    formData.project = {
      name: selectedWorkflows[workflow].inputs.accessions.display,
      desc: selectedWorkflows[workflow].inputs.accessions.display,
      type: workflow,
    }

    // set workflow inputs
    let myWorkflow = { name: workflow, input: {} }
    // set workflow input display
    let inputDisplay = { workflow: workflowList[workflow].label, input: {} }

    Object.keys(selectedWorkflows[workflow].inputs).forEach((key) => {
      myWorkflow.input[key] = selectedWorkflows[workflow].inputs[key].value
      if (selectedWorkflows[workflow].inputs[key].display) {
        inputDisplay.input[key] = selectedWorkflows[workflow].inputs[key].display
      } else {
        inputDisplay.input[selectedWorkflows[workflow].inputs[key].text] =
          selectedWorkflows[workflow].inputs[key].value
      }
    })

    // set form data
    formData.workflow = myWorkflow
    formData.inputDisplay = inputDisplay
    console.log(formData, apis.userProjects)
    // submit to server via api
    postData(apis.userProjects, formData)
      .then((data) => {
        setSubmitting(false)
        notify('success', 'Your request was submitted successfully!')
        setResetText(resetText + 1)
      })
      .catch((error) => {
        setSubmitting(false)
        alert(error)
      })
  }

  const closeMsgModal = () => {
    setOpenDialog(false)
  }

  useEffect(() => {
    setRequestSubmit(true)
    if (
      !workflow ||
      !selectedWorkflows[workflow] ||
      (selectedWorkflows[workflow] && !selectedWorkflows[workflow].validForm)
    ) {
      setRequestSubmit(false)
    }
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    setDoValidation(doValidation + 1)
  }, [workflow]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    let url = apis.userInfo
    getData(url)
      .then((data) => {
        if (data.info.allowNewRuns) {
          setDisabled(false)
        } else {
          setSysMsg(data.info.message)
          setDisabled(true)
          setOpenDialog(true)
        }
      })
      .catch((err) => {
        alert(err)
      })
  }, [])

  return (
    <div className="animated fadeIn">
      <MessageDialog
        className="modal-lg modal-danger"
        title="System Message"
        isOpen={openDialog}
        html={true}
        message={'<div><b>' + sysMsg + '</b></div>'}
        handleClickClose={closeMsgModal}
      />
      <ToastContainer />
      <LoaderDialog loading={submitting === true} text="Submitting..." />
      <Form
        onSubmit={(e) => {
          e.preventDefault()
        }}
      >
        <div className="clearfix" style={disabled ? { pointerEvents: 'none', opacity: '0.4' } : {}}>
          <h4 className="pt-3">Retrieve SRA Data</h4>
          {workflow && (
            <>
              <HtmlText
                text={
                  workflowList[workflow].info +
                  ' <a target="_blank" href=' +
                  workflowList[workflow].link +
                  ' rel="noopener noreferrer">Learn more</a>'
                }
              />
              <br></br>
            </>
          )}
          <br></br>
          {workflow === 'sra2fastq' && (
            <Sra2fastq
              reset={resetText}
              name={workflow}
              full_name={workflow}
              setParams={setWorkflowParams}
              isValid={selectedWorkflows[workflow] ? selectedWorkflows[workflow].validForm : false}
              errMessage={
                selectedWorkflows[workflow] ? selectedWorkflows[workflow].errMessage : null
              }
            />
          )}
          <br></br>
        </div>

        <div className="edge-center">
          <Button
            color="primary"
            onClick={(e) => onSubmit()}
            disabled={!workflow || !requestSubmit}
          >
            Submit
          </Button>{' '}
        </div>
        <br></br>
        <br></br>
      </Form>

      <br></br>
      <br></br>
      <SraDataTable tableType="user" title={'My SRA Data'} {...props} />
    </div>
  )
}

export default Main
