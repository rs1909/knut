#include <QDialog>
#include <QTextEdit>
#include <QScrollBar>
#include <QVBoxLayout>
#include <string>
#include <iostream>

class screenDialog : public QDialog
{
    Q_OBJECT
  public:
    screenDialog(QWidget *parent = 0) : QDialog(parent)
    {
      display = new QTextEdit;
      QVBoxLayout *layout = new QVBoxLayout;
      layout->addWidget(display);
      layout->setMargin(0);
      this->setLayout(layout);
      display->setLineWrapMode(QTextEdit::NoWrap);
      display->setReadOnly(true);
      display->setFontFamily("Courier");
      setText(" ");
//  QFontMetricsF metrics( display->document()->defaultFont() );
//  display->document()->rootFrame()->setWidth( 60*metrics.maxWidth()+2*display->frameWidth() );
//  display->document()->rootFrame()->setHeight( 20*metrics.lineSpacing()+2*display->frameWidth() );
    }
  public slots:
    void setTitle(const std::string& str)
    {
      display->setDocumentTitle(QString(str.c_str()));
    }
    void append(const std::string& str)
    {
      QTextCursor cur = display->textCursor ();
      cur.movePosition(QTextCursor::End);
      QScrollBar* bar = display->verticalScrollBar();
      const bool atEnd = (bar->maximum() == bar->value());
      cur.insertText(QString(str.c_str()));
      if (atEnd) bar->setValue(bar->maximum());
    }
    void setText(const QString& str)
    {
      display->setPlainText(str);
      display->setFontFamily("Courier");
    }
  private:
    QTextEdit *display;
};
