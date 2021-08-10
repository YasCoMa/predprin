import smtp
import os
class Communication:
	def send_success_email(task: luigi.Task) -> None:
		fromaddr = "ycfrenchgirl2@gmail.com"
        toaddr = os.environ['email_owner']
        msg = MIMEMultipart()
        msg['From'] = fromaddr
        msg['To'] = toaddr
        name_task_class = task.task_family
        msg['Subject'] = name_task_class+"ran successfully"
        body = "<p>\n"
		body += " ".join((
			name_task_class,
			"task finished successfully at",
			datetime.datetime.now().strftime("%H:%M:%S on %Y-%m-%d"),
			))
		body += "\n</p>\n"
		# Then start mining task attributes, pick some off, and dump them to a table
		body += "<table align=\"left\" border=1>\n"
		body += "\n".join([
			"<tr><td>{}</td><td>{}</td></tr>".format(key, value)
			for key, value in (
				("Task ID", task.task_id),
				("Task Class Name", name_task_class),
				("Parameter(s)", task.get_params()),
				("Requirement(s)", task.requires()),
				("Input(s)", task.input()),
				("Output(s)", task.output()),
				)
			])
		body += "\n</table>\n"
        msg.attach(MIMEText(body, 'html'))

        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.starttls()
        server.login(fromaddr, "password")
        text = msg.as_string()
        server.sendmail(fromaddr, toaddr, text)
        server.quit()

		return None
