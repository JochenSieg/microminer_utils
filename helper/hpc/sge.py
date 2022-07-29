import logging
import subprocess
import time
from pathlib import Path

logger = logging.getLogger(__name__)


class SGEJob:
    """Wraps an HPC job.

    Manages job deletion with "with"-statement.
    """

    def __init__(self, job_id: str):
        """Create new SGEJob wrapper.

        :param job_id: The jobs ID.
        """
        self.job_id = job_id

    def __enter__(self):
        """Determines what happens when with-scope is entered.

        :return: this.
        """
        # self.job_id = SGEJobRunner.submit(self.job_script)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Determines what happens when with-scope is exited.

        Clean up running jobs.

        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """
        SGEJobRunner.delete_hpc_job(self.job_id)


class SGEJobRunner:
    """Provides functionality to run and stop SGE jobs."""

    @staticmethod
    def submit_and_wait(job_script: Path) -> None:
        """Submits the given job_script to SGE and waits for completion.

        :param job_script: Path to the job script.
        :return: None
        """
        with SGEJobRunner.submit(job_script) as job:
            assert job.job_id is not None

            # wait for jobs to complete
            cmd_call = ["qstat", "-j", job.job_id]
            while True:
                time.sleep(2)
                res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
                logger.info(f"WAITING FOR HPC JOBS. QSTAT: {res.returncode}")
                if res.returncode != 0:
                    break

    @staticmethod
    def submit(job_script: Path) -> SGEJob:
        """Submit a job script to SGE.

        :param job_script: Path to the job script.
        :return: A new SGEJob instance.
        """
        cmd_call = ["qsub", str(job_script.resolve())]
        logger.info("Calling CMD: {}".format(" ".join(cmd_call)))
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        if res.returncode != 0:
            logger.error(res)
            raise ValueError("Encountered problem in job submission{}".format(res))
        hpc_jobid = None
        res_stdout = res.stdout.decode("utf-8")
        if res_stdout.startswith("Your job-array "):
            hpc_jobid = res_stdout[15:].split(".")[0]
        if hpc_jobid is None:
            raise ValueError("Unexpected return string: ", res_stdout)
        logger.info(f"Submitted job: {res_stdout}")
        return SGEJob(hpc_jobid)

    @staticmethod
    def delete_hpc_job(job_id: str) -> None:
        """Delete all SGE jobs with the provided job id.

        :param job_id: The job id of jobs to delete.
        :return: None
        """
        cmd_call = ["qdel", "-j", job_id]
        res = subprocess.run(cmd_call, stdout=subprocess.PIPE)
        # if res.returncode != 0:
        #     logger.error(res)
        #     raise ValueError('Encountered problem {}'.format(res))
