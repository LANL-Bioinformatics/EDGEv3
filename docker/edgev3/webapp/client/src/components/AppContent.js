import React, { Suspense } from 'react'
import { Navigate, Route, Routes } from 'react-router-dom'
import { CContainer, CSpinner } from '@coreui/react'
import { Row, Col } from 'reactstrap'

// routes config
import routes from 'src/routes'
import privateRoutes from 'src/private-routes'
import PrivateRoute from 'src/edge/common/PrivateRoute'

const AppContent = () => {
  return (
    <CContainer className="px-4" lg>
      <Suspense fallback={<CSpinner color="primary" />}>
        <Row className="justify-content-center">
          <Col xs="12" md="11" lg="11">
            <Routes>
              {routes.map((route, idx) => {
                return (
                  route.element && (
                    <Route
                      key={idx}
                      path={route.path}
                      exact={route.exact}
                      name={route.name}
                      element={<route.element />}
                    />
                  )
                )
              })}
              {privateRoutes.map((route, idx) => {
                return (
                  route.element && (
                    <Route
                      key={idx}
                      path={route.path}
                      exact={route.exact}
                      name={route.name}
                      element={
                        <PrivateRoute>
                          <route.element />
                        </PrivateRoute>
                      }
                    />
                  )
                )
              })}

              <Route path="*" element={<Navigate to="/home" replace />} />
            </Routes>
          </Col>
        </Row>
      </Suspense>
    </CContainer>
  )
}

export default React.memo(AppContent)
