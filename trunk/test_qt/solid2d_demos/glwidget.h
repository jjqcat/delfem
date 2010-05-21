#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

#include "delfem/camera.h"
#include "delfem/field.h"
#include "delfem/field_world.h"
#include "delfem/drawer_field.h"
#include "delfem/eqnsys_solid.h"

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();

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
private:
    Com::View::CCamera camera;

    Fem::Field::CFieldWorld world;
    Fem::Field::View::CDrawerArrayField drawer_ary;
    double cur_time;
    double dt;
    Fem::Eqn::CEqnSystem_Solid2D solid;
    unsigned int id_field_disp;
    unsigned int id_field_equiv_stress;	// Åe?Åg?ÅÒ??I?I?e
    unsigned int id_field_stress;	// Åe?Åg?ÅÒ??I?I?e
};

#endif
