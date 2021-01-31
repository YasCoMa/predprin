import luigi
from datetime import datetime
import os
class TimeTaskMixin(object):
	
	'''
	A mixin that when added to a luigi task, will print out
	the tasks execution time to standard out, when the task is
	finished
	'''
	@luigi.Task.event_handler(luigi.Event.START)
	def print_start_time(self):
		print('### START TIME ###: ' + str(datetime.now().strftime("%d/%m/%Y, %H:%M:%S")) )

	@luigi.Task.event_handler(luigi.Event.SUCCESS)
	def print_end_time(self):
		print('### END TIME ###: ' + str(datetime.now().strftime("%d/%m/%Y, %H:%M:%S")))

	@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
	def print_execution_time(self, processing_time):
		print('### PROCESSING TIME ###: ' + str(processing_time))