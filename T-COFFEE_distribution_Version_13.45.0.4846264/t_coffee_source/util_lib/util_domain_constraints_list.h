/******************************COPYRIGHT NOTICE*******************************/
/*  (c) Centro de Regulacio Genomica                                                        */
/*  and                                                                                     */
/*  Cedric Notredame                                                                        */
/*  15 Oct 2020 - 17:52.                                                                    */
/*All rights reserved.                                                                      */
/*This file is part of T-COFFEE.                                                            */
/*                                                                                          */
/*    T-COFFEE is free software; you can redistribute it and/or modify                      */
/*    it under the terms of the GNU General Public License as published by                  */
/*    the Free Software Foundation; either version 2 of the License, or                     */
/*    (at your option) any later version.                                                   */
/*                                                                                          */
/*    T-COFFEE is distributed in the hope that it will be useful,                           */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of                        */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         */
/*    GNU General Public License for more details.                                          */
/*                                                                                          */
/*    You should have received a copy of the GNU General Public License                     */
/*    along with Foobar; if not, write to the Free Software                                 */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA             */
/*...............................................                                           */
/*  If you need some more information                                                       */
/*  cedric.notredame@gmail.com                                                             */
/*...............................................                                           */
/******************************COPYRIGHT NOTICE*******************************/
Constraint_list * mask_list_with_aln (Alignment *A,int start, int len,Constraint_list *CL, int new_value);
Constraint_list* mask_list_with_aln_pair (Alignment *A,int start, int end,Constraint_list *CL,int new_value);
Constraint_list *mask_entry( Constraint_list *CL, int p, int new_value);
Constraint_list *prepare_list_and_seq4sw(Constraint_list *I, int n_seq, char **seq_name);
int ** get_undefined_list (Constraint_list *CL);
int      is_never_undefined (Constraint_list *CL,int r);
int* do_analyse_list ( Constraint_list *CL);


