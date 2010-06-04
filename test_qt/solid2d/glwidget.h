#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

#include "delfem/camera.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/eqnsys_solid.h"

class GLWidget_Solid2d : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget_Solid2d(QWidget *parent = 0);
    ~GLWidget_Solid2d();

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

private slots:
    void StepTime();

public:
    void SetNewProblem();
public:
    Com::View::CCamera camera;

    Fem::Field::CFieldWorld world;
    Fem::Field::View::CDrawerArrayField drawer_ary;
    double cur_time;
    double dt;
    Fem::Eqn::CEqnSystem_Solid2D solid;
    unsigned int id_field_disp;
    unsigned int id_field_equiv_stress;
    unsigned int id_field_stress;
};

#endif
