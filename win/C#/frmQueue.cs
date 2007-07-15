using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.Threading;
using System.Diagnostics;

namespace Handbrake
{
    public partial class frmQueue : Form
    {
        private delegate void ProgressUpdateHandler(int progressSplit);

        public frmQueue()
        {
            InitializeComponent();
        }

        private void btn_Close_Click(object sender, EventArgs e)
        {
            this.Hide();
        }

        private void btn_delete_Click(object sender, EventArgs e)
        {
            list_queue.Items.Remove(list_queue.SelectedItem);
        }

        private void btn_up_Click(object sender, EventArgs e)
        {
            int count = list_queue.Items.Count;
            int itemToMove = list_queue.SelectedIndex;
            int previousItemint = 0; 
            String previousItem = "";

             if (itemToMove > 0){
                previousItemint = itemToMove - 1;
                previousItem = list_queue.Items[previousItemint].ToString();
                list_queue.Items[previousItemint] = list_queue.Items[itemToMove];
                list_queue.Items[itemToMove] = previousItem;
                list_queue.SelectedIndex = list_queue.SelectedIndex - 1;
            }
        }

        private void btn_down_Click(object sender, EventArgs e)
        {
            int count = list_queue.Items.Count;
            int itemToMove = list_queue.SelectedIndex;
            int itemAfterInt = 0; 
            String itemAfter = "";

            if (itemToMove < (count - 1))
            {
                itemAfterInt = itemToMove + 1;
                itemAfter = list_queue.Items[itemAfterInt].ToString();
                list_queue.Items[itemAfterInt] = list_queue.Items[itemToMove];
                list_queue.Items[itemToMove] = itemAfter;
                list_queue.SelectedIndex = list_queue.SelectedIndex + 1;
            }
        }

        private void btn_q_encoder_Click(object sender, EventArgs e)
        {
            progressBar.Value = 0;
            lbl_progressValue.Text = "0 %";
            progressBar.Update();
            ThreadPool.QueueUserWorkItem(startProc);
        }

        private void startProc(object state)
        {
            int listSize = list_queue.Items.Count;
            listSize--;
            int counter = 0;
            int progressSplit = listSize / 100;

            while (counter <= listSize)
            {
                String query = list_queue.Items[0].ToString();

                Process hbProc = new Process();
                hbProc.StartInfo.FileName = "hbcli.exe";
                hbProc.StartInfo.Arguments = query;
                hbProc.StartInfo.UseShellExecute = false;
                hbProc.Start();

                // Set the process Priority
                string priority = Properties.Settings.Default.processPriority;
                switch (priority)
                {
                    case "Realtime":
                        hbProc.PriorityClass = ProcessPriorityClass.RealTime;
                        break;
                    case "High":
                        hbProc.PriorityClass = ProcessPriorityClass.High;
                        break;
                    case "Above Normal":
                        hbProc.PriorityClass = ProcessPriorityClass.AboveNormal;
                        break;
                    case "Normal":
                        hbProc.PriorityClass = ProcessPriorityClass.Normal;
                        break;
                    case "Low":
                        hbProc.PriorityClass = ProcessPriorityClass.Idle;
                        break;
                    default:
                        hbProc.PriorityClass = ProcessPriorityClass.BelowNormal;
                        break;
                }

                hbProc.WaitForExit();
                hbProc.Close();
                counter++;

                updateUIElements(progressSplit);
            }
        }

        private void updateUIElements(int progressSplit)
        {
            if (this.InvokeRequired)
            {
                this.BeginInvoke(new ProgressUpdateHandler(updateUIElements), new object[] { progressSplit });
                return;
            }

            this.list_queue.Items.Remove(0);
            progressBar.Value = progressBar.Value + progressSplit;
            progressBar.Update();
        }
    }
}