using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    internal class RentalType
    {
        [Key]
        public int Id { get; set; }

        [ForeignKey("RentalId")]
        public int RentalId { get; set; }


        [Required]
        public string? Description { get; set; } //description of type
    }
}
