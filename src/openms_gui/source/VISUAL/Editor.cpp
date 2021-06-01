
#include <OpenMS/VISUAL/Editor.h>
#include <ui_Editor.h>

#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <QMessageBox>

namespace OpenMS {
  Editor::Editor(QWidget* parent) :
      QWidget(parent), ui(new Ui::Editor)
  {
    ui->setupUi(this);
    QObject::connect(ui->New, SIGNAL(clicked()), this, SLOT(on_New_clicked()));
    QObject::connect(ui->Save, SIGNAL(clicked()), this, SLOT(on_Save_clicked()));
    QObject::connect(ui->SaveAs, SIGNAL(clicked()), this, SLOT(on_Save_As_clicked()));
    QObject::connect(ui->Copy, SIGNAL(clicked()), this, SLOT(on_Copy_clicked()));
    QObject::connect(ui->Paste, SIGNAL(clicked()), this, SLOT(on_Paste_clicked()));
    QObject::connect(ui->Undo, SIGNAL(clicked()), this, SLOT(on_Undo_clicked()));
    QObject::connect(ui->Redo, SIGNAL(clicked()), this, SLOT(on_Redo_clicked()));
    QObject::connect(ui->Cut, SIGNAL(clicked()), this, SLOT(on_Cut_clicked()));
    QObject::connect(ui->Open, SIGNAL(clicked()), this, SLOT(on_Open_clicked()));
  }

  Editor::~Editor()
  {
    delete ui;
  }

  void Editor::on_New_clicked()
  {
    file_path_ = "";
    ui->textEdit->setText("");
  }

  void Editor::on_Open_clicked()
  {
    QString file_name = QFileDialog::getOpenFileName(this, "Open the file");
    QFile file(file_name);
    file_path_ = file_name;
    if (!file.open(QFile::ReadOnly | QFile::Text))
    {
      QMessageBox::warning(this, "..", "file not open");
      return;
    }
    QTextStream in(&file);
    QString text = in.readAll();
    ui->textEdit->setText(text);
    file.close();
  }

  void Editor::on_Save_clicked()
  {
    QFile file(file_path_);
    if (!file.open(QFile::ReadOnly | QFile::Text))
    {
      QMessageBox::warning(this, "..", "file not open");
      return;
    }
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.flush();
    file.close();
  }

  void Editor::on_Save_As_clicked()
  {
    QString file_name = QFileDialog::getSaveFileName(this, "Open the file");
    QFile file(file_name);
    file_path_ = file_name;
    if (!file.open(QFile::WriteOnly | QFile::Text))
    {
      QMessageBox::warning(this, "..", "file not open");
      return;
    }
    QTextStream out(&file);
    QString text = ui->textEdit->toPlainText();
    out << text;
    file.flush();
    file.close();
  }

  void Editor::on_Copy_clicked()
  {
    ui->textEdit->copy();
  }

  void Editor::on_Paste_clicked()
  {
    ui->textEdit->paste();
  }

  void Editor::on_Cut_clicked()
  {
    ui->textEdit->cut();
  }

  void Editor::on_Undo_clicked()
  {
    ui->textEdit->undo();
  }

  void Editor::on_Redo_clicked()
  {
    ui->textEdit->redo();
  }
}
