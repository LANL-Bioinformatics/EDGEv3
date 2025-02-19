import React, { useState, useEffect } from 'react'
import { Button, Form } from 'reactstrap'
import { useNavigate } from 'react-router-dom'
import { ToastContainer } from 'react-toastify'
import 'react-toastify/dist/ReactToastify.css'

import { workflowList } from 'src/util'
import { postData, getData, notify, apis } from '../../common/util'
import { LoaderDialog, MessageDialog } from '../../common/Dialogs'
import MySelect from '../../common/MySelect'
import { Project } from '../../project/forms/Project'
import { HtmlText } from '../../common/HtmlText'
import { RunFaQCs } from './forms/RunFaQCs'
import { workflowOptions } from './defaults'

const Main = (props) => {
  const navigate = useNavigate()
  const [submitting, setSubmitting] = useState(false)
  const [requestSubmit, setRequestSubmit] = useState(false)
  const [projectParams, setProjectParams] = useState()
  const [selectedWorkflows, setSelectedWorkflows] = useState({})
  const [doValidation, setDoValidation] = useState(0)
  const [workflow, setWorkflow] = useState(workflowOptions[0].value)
  const [openDialog, setOpenDialog] = useState(false)
  const [disabled, setDisabled] = useState(false)
  const [sysMsg, setSysMsg] = useState()

  //callback function for child component
  const setProject = (params) => {
    //console.log('main project:', params)
    setProjectParams(params)
    setDoValidation(doValidation + 1)
  }
  const setWorkflowParams = (params, workflowName) => {
    //console.log(workflowName, params)
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
      name: projectParams.projectName,
      desc: projectParams.projectDesc,
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
    // set pairedFile flag
    myWorkflow.input['pairedFile'] = !myWorkflow.input['interleaved']

    //update input for nextflow
    if (!myWorkflow.input['artifactFile']) {
      myWorkflow.input['artifactFile'] = '${projectDir}/nf_assets/NO_FILE3'
    }
    // set form data
    formData.workflow = myWorkflow
    formData.inputDisplay = inputDisplay

    // files used for caculating total input size on server side
    formData.files = selectedWorkflows[workflow].files

    // submit to server via api
    postData(apis.userProjects, formData)
      .then((data) => {
        setSubmitting(false)
        notify('success', 'Your workflow request was submitted successfully!', 2000)
        setTimeout(() => navigate('/user/projects'), 2000)
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

    if (projectParams && !projectParams.validForm) {
      setRequestSubmit(false)
    }

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
  }, [props])

  return (
    <div
      className="animated fadeIn"
      style={disabled ? { pointerEvents: 'none', opacity: '0.4' } : {}}
    >
      <MessageDialog
        className="modal-lg modal-danger"
        title="System Message"
        isOpen={openDialog}
        html={true}
        message={'<div><b>' + sysMsg + '</b></div>'}
        handleClickClose={closeMsgModal}
      />
      <span className="pt-3 text-muted edge-text-size-small">
        Metagenomics | Run Single Workflow{' '}
      </span>
      <ToastContainer />
      <LoaderDialog loading={submitting === true} text="Submitting..." />
      <Form
        onSubmit={(e) => {
          e.preventDefault()
        }}
      >
        <div className="clearfix">
          <h4 className="pt-3">Run a Single Workflow</h4>
          <hr />
          <Project setParams={setProject} />

          <br></br>
          <b>Workflow</b>
          <MySelect
            //defaultValue={workflowOptions[0]}
            options={workflowOptions}
            value={workflowOptions[0]}
            onChange={(e) => {
              if (e) {
                setWorkflow(e.value)
              } else {
                setWorkflow()
              }
            }}
            placeholder="Select a Workflow..."
            isClearable={true}
          />
          {workflow && (
            <div className="pt-3 text-muted edge-text-size-small">
              <HtmlText
                text={
                  workflowList[workflow].info +
                  ' <a target="_blank" href=' +
                  workflowList[workflow].link +
                  ' rel="noopener noreferrer">Learn more</a>'
                }
              />
              <br></br>
            </div>
          )}
          <br></br>
          {workflow === 'runFaQCs' && (
            <RunFaQCs
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
    </div>
  )
}

export default Main
