//todo, actually make these into a good unit test for the partial_conf class
  const int num_modes = 300;
  const int num_bits = 4;
  const int giant_words = (num_bits*num_modes)/64 +1;
  partial_config<num_modes,num_bits,giant_words> c;
  //c.print_members();
  c.set_mode(0,3);
  c.set_mode(298,3);
  c.set_mode(299, 5);
  c.print_n_words();
  c.decrement_mode(0);
  //c.increment_mode(298);
  //c.print_n_words();
  for(int i = 0; i < num_modes; i++){
    for(int j = 0; j < 16; j++){
      c.set_mode(i,j);
      cout << "conf mode "<< i <<" : " << c.get_mode(i) << " ?= " << j << endl;
      assert( c.get_mode(i) == j);
    }
  }
