import React, { useEffect, useRef } from 'react'
import { NavLink, useNavigate } from 'react-router-dom'
import { useSelector, useDispatch } from 'react-redux'
import { setSidebar } from '../redux/reducers/coreuiSlice'
import {
  CContainer,
  CDropdown,
  CDropdownItem,
  CDropdownMenu,
  CDropdownToggle,
  CHeader,
  CHeaderNav,
  CHeaderToggler,
  CNavLink,
  CNavItem,
  useColorModes,
  CHeaderBrand,
} from '@coreui/react'
import CIcon from '@coreui/icons-react'
import {
  cilBell,
  cilContrast,
  cilEnvelopeOpen,
  cilList,
  cilMenu,
  cilMoon,
  cilSun,
} from '@coreui/icons'

import { AppHeaderDropdown } from './header/index'
import logo from 'src/assets/brand/logo.png'
import { logout } from 'src/redux/reducers/edge/userSlice'

const AppHeader = () => {
  const headerRef = useRef()
  const { colorMode, setColorMode } = useColorModes('coreui-free-react-admin-template-theme')

  const dispatch = useDispatch()
  const sidebarShow = useSelector((state) => state.coreui.sidebarShow)
  const navigate = useNavigate()
  const user = useSelector((state) => state.user)

  const signOut = (e) => {
    dispatch(logout())
    navigate('/login')
  }

  useEffect(() => {
    document.addEventListener('scroll', () => {
      headerRef.current &&
        headerRef.current.classList.toggle('shadow-sm', document.documentElement.scrollTop > 0)
    })
  }, [])

  return (
    <CHeader position="sticky" className="mb-4 p-0" ref={headerRef}>
      <CContainer className="border-bottom px-4" fluid>
        {!sidebarShow && (
          <CHeaderBrand className="sidebar-brand-narrow" to="/">
            <img alt="logo" style={{ width: 150, height: 50 }} src={logo} />
          </CHeaderBrand>
        )}
        <CHeaderToggler
          onClick={() => dispatch(setSidebar(!sidebarShow))}
          style={{ marginInlineStart: '-14px' }}
        >
          <CIcon icon={cilMenu} size="lg" />
        </CHeaderToggler>

        {user.isAuthenticated ? (
          <>
            <CHeaderNav className="d-none d-md-flex me-auto">
              <CNavItem>
                <CNavLink to="/user/projects" as={NavLink}>
                  My Projects
                </CNavLink>
              </CNavItem>
              <CNavItem>
                <CNavLink to="/user/uploads" as={NavLink}>
                  My Uploads
                </CNavLink>
              </CNavItem>
              {/*  <CNavItem>
                <CNavLink to="/user/sradata" as={NavLink}>
                  My SRA Data
                </CNavLink>
              </CNavItem> */}
              <CNavItem>
                <CNavLink to="/user/jobqueue" as={NavLink}>
                  Job Queue
                </CNavLink>
              </CNavItem>
            </CHeaderNav>
            <CHeaderNav className="ms-3">
              <AppHeaderDropdown user={user} logout={(e) => signOut(e)} />
            </CHeaderNav>
          </>
        ) : (
          <>
            <CHeaderNav className="d-none d-md-flex"></CHeaderNav>
            <CHeaderNav className="ms-auto">
              <CNavLink to="/login" as={NavLink}>
                Login
              </CNavLink>
              <CNavLink to="/register" as={NavLink}>
                Sign up
              </CNavLink>
            </CHeaderNav>
          </>
        )}
        {/* <CHeaderNav>
          <li className="nav-item py-1">
            <div className="vr h-100 mx-2 text-body text-opacity-75"></div>
          </li>
          <CDropdown variant="nav-item" placement="bottom-end">
            <CDropdownToggle caret={false}>
              {colorMode === 'dark' ? (
                <CIcon icon={cilMoon} size="lg" />
              ) : colorMode === 'auto' ? (
                <CIcon icon={cilContrast} size="lg" />
              ) : (
                <CIcon icon={cilSun} size="lg" />
              )}
            </CDropdownToggle>
            <CDropdownMenu>
              <CDropdownItem
                active={colorMode === 'light'}
                className="d-flex align-items-center"
                as="button"
                type="button"
                onClick={() => setColorMode('light')}
              >
                <CIcon className="me-2" icon={cilSun} size="lg" /> Light
              </CDropdownItem>
              <CDropdownItem
                active={colorMode === 'dark'}
                className="d-flex align-items-center"
                as="button"
                type="button"
                onClick={() => setColorMode('dark')}
              >
                <CIcon className="me-2" icon={cilMoon} size="lg" /> Dark
              </CDropdownItem>
              <CDropdownItem
                active={colorMode === 'auto'}
                className="d-flex align-items-center"
                as="button"
                type="button"
                onClick={() => setColorMode('auto')}
              >
                <CIcon className="me-2" icon={cilContrast} size="lg" /> Auto
              </CDropdownItem>
            </CDropdownMenu>
          </CDropdown>
          <li className="nav-item py-1">
            <div className="vr h-100 mx-2 text-body text-opacity-75"></div>
          </li>
        </CHeaderNav> */}
      </CContainer>
    </CHeader>
  )
}

export default AppHeader
