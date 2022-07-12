from pathlib import Path
import time
import subprocess
import logging

logger = logging.getLogger(__name__)


class SGEJob:
    """
    Wraps a HPC job. Manages job deletion with "with"-statement.
    """

    def __init__(self, job_id: str):
        self.job_id = job_id

    def __enter__(self):
        # self.job_id = SGEJobRunner.submit(self.job_script)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        SGEJobRunner.delete_hpc_job(self.job_id)


class SGEJobRunner:

    @staticmethod
    def submit_and_wait(job_script: Path):
        with SGEJobRunner.submit(job_script) as job:
            assert job.job_id is not None

            # wait for jobs to complete
            cmd_call = ['qstat', '-j', job.job_id]
            while True:
                time.sleep(2)
                res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
                logger.info(f'WAITING FOR HPC JOBS. QSTAT: {res.returncode}')
                if res.returncode != 0:
                    break

                #     logger.error(res)
                #     raise ValueError('Encountered problem {}'.format(res))
                # lines = res.stdout.decode().split('\n')
                # print('QSTAT', lines)
                # if not any(l.strip().startswith(hpc_jobid) for l in lines):
                #     break

    @staticmethod
    def submit(job_script: Path):
        cmd_call = ['qsub', str(job_script.resolve())]
        logger.info('Calling CMD: {}'.format(' '.join(cmd_call)))
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        if res.returncode != 0:
            logger.error(res)
            raise ValueError('Encountered problem in job submission{}'.format(res))
        hpc_jobid = None
        res_stdout = res.stdout.decode('utf-8')
        if res_stdout.startswith('Your job-array '):
            hpc_jobid = res_stdout[15:].split('.')[0]
        if hpc_jobid is None:
            raise ValueError('Unexpected return string: ', res_stdout)
        logger.info(f'Submitted job: {res_stdout}')
        return SGEJob(hpc_jobid)

    @staticmethod
    def delete_hpc_job(job_id: str):
        cmd_call = ['qdel', '-j', job_id]
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        # if res.returncode != 0:
        #     logger.error(res)
        #     raise ValueError('Encountered problem {}'.format(res))

