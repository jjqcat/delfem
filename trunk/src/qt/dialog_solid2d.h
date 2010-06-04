#ifndef DIALOG_SOLID2D_H
#define DIALOG_SOLID2D_H

#include <QDialog>

namespace Ui {
    class Dialog_Solid2D;
}

class Dialog_Solid2D : public QDialog {
    Q_OBJECT
public:
    Dialog_Solid2D(QWidget *parent = 0);
    ~Dialog_Solid2D();

protected:
    void changeEvent(QEvent *e);

private:
    Ui::Dialog_Solid2D *ui;
};

#endif // DIALOG_SOLID2D_H
