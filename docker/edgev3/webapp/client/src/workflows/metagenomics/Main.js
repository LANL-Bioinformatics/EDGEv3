import React, { useState, useEffect } from 'react'
import { Button, Form } from 'reactstrap'
import { useNavigate } from 'react-router-dom'
import { ToastContainer } from 'react-toastify'
import 'react-toastify/dist/ReactToastify.css'

import { workflowList } from 'src/util'
import { postData, getData, notify, apis, isValidFileInput } from 'src/edge/common/util'
import { LoaderDialog, MessageDialog } from 'src/edge/common/Dialogs'
import MySelect from 'src/edge/common/MySelect'
import { Project } from 'src/edge/project/forms/Project'
import { HtmlText } from 'src/edge/common/HtmlText'
import { InputRawReads } from './forms/InputRawReads'
import { RunFaQCs } from './forms/RunFaQCs'
import { Assembly } from './forms/Assembly'
import { Annotation } from './forms/Annotation'
import { Binning } from './forms/Binning'
import { workflowOptions, workflows } from './defaults'

const Main = (props) => {
  const navigate = useNavigate()
  const [submitting, setSubmitting] = useState(false)
  const [requestSubmit, setRequestSubmit] = useState(false)
  const [projectParams, setProjectParams] = useState()
  const [rawDataParams, setRawDataParams] = useState()
  const [selectedWorkflows, setSelectedWorkflows] = useState({})
  const [doValidation, setDoValidation] = useState(0)
  const [workflow, setWorkflow] = useState(workflowOptions[0].value)
  const [openDialog, setOpenDialog] = useState(false)
  const [disabled, setDisabled] = useState(false)
  const [sysMsg, setSysMsg] = useState()
  const [allExpand, setAllExpand] = useState(0)
  const [allClosed, setAllClosed] = useState(0)
  //disable the expand | close
  const disableExpandClose = false

  //callback function for child component
  const setProject = (params) => {
    //console.log('main project:', params)
    setProjectParams(params)
    setDoValidation(doValidation + 1)
  }
  //callback function for child component
  const setRawData = (params) => {
    //console.log('rawData:', params)
    setRawDataParams(params)
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
    if (workflow === 'annotation' || workflow === 'binning') {
      formData.rawReads = {
        source: rawDataParams.inputs.source.value,
        inputFasta: rawDataParams.inputs.inputFiles.value[0],
      }
    } else {
      formData.rawReads = {
        source: rawDataParams.inputs.source.value,
        seqPlatform: rawDataParams.inputs.seqPlatform.value,
        paired: rawDataParams.inputs.paired.value,
        inputFiles: rawDataParams.inputs.inputFiles.value,
      }
    }

    // set workflow inputs
    let myWorkflow = { name: workflow, input: {} }
    // set workflow input display
    let inputDisplay = { 'Raw Reads': {} }
    inputDisplay[workflowList[workflow].label] = {}
    if (workflow === 'annotation' || workflow === 'binning') {
      inputDisplay['Raw Reads'][rawDataParams.inputs['source'].text] =
        rawDataParams.inputs['source'].display
      inputDisplay['Raw Reads']['Contig/Fasta File'] = rawDataParams.inputs['inputFiles'].display[0]
    } else {
      Object.keys(rawDataParams.inputs).forEach((key) => {
        if (rawDataParams.inputs[key].display) {
          inputDisplay['Raw Reads'][rawDataParams.inputs[key].text] =
            rawDataParams.inputs[key].display
        } else {
          inputDisplay['Raw Reads'][rawDataParams.inputs[key].text] =
            rawDataParams.inputs[key].value
        }
      })
    }

    //add selected assembler inputs to main inputs
    if (workflow === 'assembly') {
      // eslint-disable-next-line prettier/prettier
      selectedWorkflows[workflow].inputs = {
        ...selectedWorkflows[workflow].inputs,
        // eslint-disable-next-line prettier/prettier
        ...selectedWorkflows[workflow].assemblerInputs[selectedWorkflows[workflow].inputs['assembler'].value]
      }
    }
    //add selected annotateProgram inputs to main inputs
    if (workflow === 'annotation') {
      // eslint-disable-next-line prettier/prettier
      selectedWorkflows[workflow].inputs = {
        ...selectedWorkflows[workflow].inputs,
        // eslint-disable-next-line prettier/prettier
        ...selectedWorkflows[workflow].annotateProgramInputs[selectedWorkflows[workflow].inputs['annotateProgram'].value]
      }
    }

    Object.keys(selectedWorkflows[workflow].inputs).forEach((key) => {
      myWorkflow.input[key] = selectedWorkflows[workflow].inputs[key].value
      if (selectedWorkflows[workflow].inputs[key].display) {
        inputDisplay[[workflowList[workflow].label]][selectedWorkflows[workflow].inputs[key].text] =
          selectedWorkflows[workflow].inputs[key].display
      } else {
        inputDisplay[[workflowList[workflow].label]][selectedWorkflows[workflow].inputs[key].text] =
          selectedWorkflows[workflow].inputs[key].value
      }
    })

    // set form data
    formData.workflow = myWorkflow
    formData.inputDisplay = inputDisplay

    // files used for caculating total input size on server side
    formData.files = [...rawDataParams.files, ...selectedWorkflows[workflow].files]

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
    if (rawDataParams && !rawDataParams.validForm) {
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
              setAllExpand(0)
              setAllClosed(0)
              if (e) {
                setWorkflow(e.value)
              } else {
                setWorkflow()
              }
            }}
            placeholder="Select a Workflow..."
            isClearable={true}
          />
          {workflow && workflowList[workflow].info && (
            <div className="pt-3 text-muted edge-text-size-small">
              <HtmlText text={workflowList[workflow].info} />
              <br></br>
            </div>
          )}
          <br></br>
          {!disableExpandClose && (
            <>
              <div className="float-end edge-text-size-small">
                <Button
                  style={{ fontSize: 12, paddingBottom: '5px' }}
                  size="sm"
                  className="btn-pill"
                  color="ghost-primary"
                  onClick={() => setAllExpand(allExpand + 1)}
                >
                  expand
                </Button>
                &nbsp; | &nbsp;
                <Button
                  style={{ fontSize: 12, paddingBottom: '5px' }}
                  size="sm"
                  className="btn-pill"
                  color="ghost-primary"
                  onClick={() => setAllClosed(allClosed + 1)}
                >
                  close
                </Button>
                &nbsp; all sections &nbsp;
              </div>
              <br></br>
              <br></br>
            </>
          )}
          {workflow === 'runFaQCs' && (
            <>
              <InputRawReads
                setParams={setRawData}
                //sourceOptions={workflows[workflow]['rawReadsInput'].sourceOptions}
                isValidFileInput={isValidFileInput}
                source={workflows[workflow]['rawReadsInput'].source}
                sourceDisplay={workflows[workflow]['rawReadsInput'].text}
                text={workflows[workflow]['rawReadsInput'].text}
                tooltip={workflows[workflow]['rawReadsInput'].tooltip}
                enableInput={workflows[workflow]['rawReadsInput'].enableInput}
                placeholder={workflows[workflow]['rawReadsInput'].placeholder}
                dataSources={workflows[workflow]['rawReadsInput'].dataSources}
                fileTypes={workflows[workflow]['rawReadsInput'].fileTypes}
                projectTypes={workflows[workflow]['rawReadsInput'].projectTypes}
                projectScope={workflows[workflow]['rawReadsInput'].projectScope}
                viewFile={workflows[workflow]['rawReadsInput'].viewFile}
                isOptional={workflows[workflow]['rawReadsInput'].isOptional}
                cleanupInput={workflows[workflow]['rawReadsInput'].cleanupInput}
                maxInput={workflows[workflow]['rawReadsInput'].maxInput}
                isValid={rawDataParams ? rawDataParams.validForm : false}
                errMessage={rawDataParams ? rawDataParams.errMessage : null}
                allExpand={allExpand}
                allClosed={allClosed}
              />
              <RunFaQCs
                name={workflow}
                full_name={workflow}
                setParams={setWorkflowParams}
                isValid={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].validForm : false
                }
                errMessage={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].errMessage : null
                }
                allExpand={allExpand}
                allClosed={allClosed}
              />
            </>
          )}
          {workflow === 'assembly' && (
            <>
              <InputRawReads
                setParams={setRawData}
                isValidFileInput={isValidFileInput}
                source={workflows[workflow]['rawReadsInput'].source}
                sourceDisplay={workflows[workflow]['rawReadsInput'].text}
                text={workflows[workflow]['rawReadsInput'].text}
                tooltip={workflows[workflow]['rawReadsInput'].tooltip}
                enableInput={workflows[workflow]['rawReadsInput'].enableInput}
                placeholder={workflows[workflow]['rawReadsInput'].placeholder}
                dataSources={workflows[workflow]['rawReadsInput'].dataSources}
                fileTypes={workflows[workflow]['rawReadsInput'].fileTypes}
                projectTypes={workflows[workflow]['rawReadsInput'].projectTypes}
                projectScope={workflows[workflow]['rawReadsInput'].projectScope}
                viewFile={workflows[workflow]['rawReadsInput'].viewFile}
                isOptional={workflows[workflow]['rawReadsInput'].isOptional}
                cleanupInput={workflows[workflow]['rawReadsInput'].cleanupInput}
                maxInput={workflows[workflow]['rawReadsInput'].maxInput}
                isValid={rawDataParams ? rawDataParams.validForm : false}
                errMessage={rawDataParams ? rawDataParams.errMessage : null}
                allExpand={allExpand}
                allClosed={allClosed}
              />
              <Assembly
                name={workflow}
                full_name={workflow}
                setParams={setWorkflowParams}
                isValid={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].validForm : false
                }
                errMessage={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].errMessage : null
                }
                seqPlatform={rawDataParams.inputs.seqPlatform.value}
                allExpand={allExpand}
                allClosed={allClosed}
              />
            </>
          )}
          {workflow === 'annotation' && (
            <>
              <InputRawReads
                setParams={setRawData}
                isValidFileInput={isValidFileInput}
                source={workflows[workflow]['rawReadsInput'].source}
                sourceDisplay={workflows[workflow]['rawReadsInput'].text}
                text={workflows[workflow]['rawReadsInput'].text}
                tooltip={workflows[workflow]['rawReadsInput'].tooltip}
                enableInput={workflows[workflow]['rawReadsInput'].enableInput}
                placeholder={workflows[workflow]['rawReadsInput'].placeholder}
                dataSources={workflows[workflow]['rawReadsInput'].dataSources}
                fileTypes={workflows[workflow]['rawReadsInput'].fileTypes}
                projectTypes={workflows[workflow]['rawReadsInput'].projectTypes}
                projectScope={workflows[workflow]['rawReadsInput'].projectScope}
                viewFile={workflows[workflow]['rawReadsInput'].viewFile}
                isOptional={workflows[workflow]['rawReadsInput'].isOptional}
                cleanupInput={workflows[workflow]['rawReadsInput'].cleanupInput}
                maxInput={workflows[workflow]['rawReadsInput'].maxInput}
                isValid={rawDataParams ? rawDataParams.validForm : false}
                errMessage={rawDataParams ? rawDataParams.errMessage : null}
                allExpand={allExpand}
                allClosed={allClosed}
              />
              <Annotation
                name={workflow}
                full_name={workflow}
                setParams={setWorkflowParams}
                isValid={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].validForm : false
                }
                errMessage={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].errMessage : null
                }
                allExpand={allExpand}
                allClosed={allClosed}
              />
            </>
          )}
          {workflow === 'binning' && (
            <>
              <InputRawReads
                setParams={setRawData}
                isValidFileInput={isValidFileInput}
                source={workflows[workflow]['rawReadsInput'].source}
                sourceDisplay={workflows[workflow]['rawReadsInput'].text}
                text={workflows[workflow]['rawReadsInput'].text}
                tooltip={workflows[workflow]['rawReadsInput'].tooltip}
                enableInput={workflows[workflow]['rawReadsInput'].enableInput}
                placeholder={workflows[workflow]['rawReadsInput'].placeholder}
                dataSources={workflows[workflow]['rawReadsInput'].dataSources}
                fileTypes={workflows[workflow]['rawReadsInput'].fileTypes}
                projectTypes={workflows[workflow]['rawReadsInput'].projectTypes}
                projectScope={workflows[workflow]['rawReadsInput'].projectScope}
                viewFile={workflows[workflow]['rawReadsInput'].viewFile}
                isOptional={workflows[workflow]['rawReadsInput'].isOptional}
                cleanupInput={workflows[workflow]['rawReadsInput'].cleanupInput}
                maxInput={workflows[workflow]['rawReadsInput'].maxInput}
                isValid={rawDataParams ? rawDataParams.validForm : false}
                errMessage={rawDataParams ? rawDataParams.errMessage : null}
                allExpand={allExpand}
                allClosed={allClosed}
              />
              <Binning
                name={workflow}
                full_name={workflow}
                setParams={setWorkflowParams}
                isValid={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].validForm : false
                }
                errMessage={
                  selectedWorkflows[workflow] ? selectedWorkflows[workflow].errMessage : null
                }
                allExpand={allExpand}
                allClosed={allClosed}
              />
            </>
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
