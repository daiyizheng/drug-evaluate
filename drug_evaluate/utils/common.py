# -*- encoding: utf-8 -*-
'''
Filename         :common.py
Description      :
Time             :2023/07/19 16:48:41
Author           :daiyizheng
Email            :387942239@qq.com
Version          :1.0
'''
from __future__ import absolute_import, print_function, annotations
import logging

from rdkit import rdBase

logger = logging.getLogger(__name__)

def disable_rdkit_log():
    rdBase.DisableLog('rdApp.*')


def enable_rdkit_log():
    rdBase.EnableLog('rdApp.*')